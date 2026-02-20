#ifndef SZ3_SPERR_CORE_IMPL_HPP
#define SZ3_SPERR_CORE_IMPL_HPP

// Header-only SPERR core implementation for SZ3 integration.
// Generated from upstream SPERR source files and inlined for ODR safety.

// ---- BEGIN sperr_helper.cpp ----
#include "sperr_helper.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

inline auto SZ3::SPERR::aligned_malloc(size_t alignment, size_t size) -> void*
{
#ifdef _WIN32
  return _aligned_malloc(size, alignment);
#else
  return std::aligned_alloc(alignment, size);
#endif
}

inline void SZ3::SPERR::aligned_free(void* p)
{
#ifdef _WIN32
  _aligned_free(p);
#else
  std::free(p);
#endif
}

inline auto SZ3::SPERR::num_of_xforms(size_t len) -> size_t
{
  assert(len > 0);
  // I decide 9 is the minimal length to do one level of xform.
  // I also decide that no matter what the input size is,
  // six (6) is the maxinum number of transforms to do.
  //
  size_t num = 0;
  while (len >= 9) {
    ++num;
    len -= len / 2;
  }
  return std::min(num, size_t{6});
}

inline auto SZ3::SPERR::can_use_dyadic(dims_type dims) -> std::optional<size_t>
{
  // In case of 2D or 1D `dims`, return empty right away.
  if (dims[2] < 2 || dims[1] < 2)
    return {};

  auto xy = SZ3::SPERR::num_of_xforms(std::min(dims[0], dims[1]));
  auto z = SZ3::SPERR::num_of_xforms(dims[2]);

  // Note: if some dimensions can do 5 levels of transforms and some can do 6, we use
  //       dyanic scheme and do 5 levels on all of them. I.e., the benefit of dyanic
  //       transforms exceeds one extra level of transform.
  //
  if ((xy == z) || (xy >= 5 && z >= 5))
    return std::min(xy, z);
  else
    return {};
}

inline auto SZ3::SPERR::coarsened_resolutions(dims_type full_dims) -> std::vector<dims_type>
{
  auto resolutions = std::vector<dims_type>();

  if (full_dims[2] > 1) {  // 3D.
    const auto dyadic = SZ3::SPERR::can_use_dyadic(full_dims);
    if (dyadic) {
      resolutions.reserve(*dyadic);
      for (size_t lev = *dyadic; lev > 0; lev--) {
        auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(full_dims[0], lev);
        auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(full_dims[1], lev);
        auto [z, zd] = SZ3::SPERR::calc_approx_detail_len(full_dims[2], lev);
        resolutions.push_back({x, y, z});
      }
    }
  }
  else {  // 2D. Assume that there's no 1D use case that requires multi-resolution.
    size_t xy = SZ3::SPERR::num_of_xforms(std::min(full_dims[0], full_dims[1]));
    resolutions.reserve(xy);
    for (size_t lev = xy; lev > 0; lev--) {
      auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(full_dims[0], lev);
      auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(full_dims[1], lev);
      resolutions.push_back({x, y, 1});
    }
  }

  return resolutions;
}

inline auto SZ3::SPERR::coarsened_resolutions(dims_type vdim, dims_type cdim) -> std::vector<dims_type>
{
  auto resolutions = std::vector<dims_type>();

  // Test if the volume dimension is divisible by the chunk dimension.
  bool divisible = true;
  for (size_t i = 0; i < 3; i++)
    if (vdim[i] % cdim[i] != 0)
      divisible = false;

  if (divisible) {
    auto nx = vdim[0] / cdim[0];
    auto ny = vdim[1] / cdim[1];
    auto nz = vdim[2] / cdim[2];

    resolutions = SZ3::SPERR::coarsened_resolutions(cdim);
    for (auto& resolution : resolutions) {
      resolution[0] *= nx;
      resolution[1] *= ny;
      resolution[2] *= nz;
    }
  }

  return resolutions;
}

inline auto SZ3::SPERR::num_of_partitions(size_t len) -> size_t
{
  size_t num_of_parts = 0;  // Num. of partitions we can do
  while (len > 1) {
    num_of_parts++;
    len -= len / 2;
  }

  return num_of_parts;
}

inline auto SZ3::SPERR::calc_approx_detail_len(size_t orig_len, size_t lev) -> std::array<size_t, 2>
{
  size_t low_len = orig_len;
  size_t high_len = 0;
  for (size_t i = 0; i < lev; i++) {
    high_len = low_len / 2;
    low_len -= high_len;
  }

  return {low_len, high_len};
}

// Good solution to deal with bools and unsigned chars
// https://stackoverflow.com/questions/8461126/how-to-create-a-byte-out-of-8-bool-values-and-vice-versa
inline auto SZ3::SPERR::pack_booleans(vec8_type& dest, const std::vector<bool>& src, size_t offset) -> RTNType
{
  // `src` has to have a size of multiples of 8.
  if (src.size() % 8 != 0)
    return RTNType::WrongLength;

  // `dest` should have enough space, as the API specifies.
  assert(dest.size() >= offset + src.size() / 8);

  // How many bits to process at a time.
  constexpr uint64_t bit_stride = 2048;
  constexpr uint64_t byte_stride = bit_stride / 8;
  constexpr uint64_t magic = 0x8040201008040201;

  // It turns out C++ doesn't specify bool to be 1 byte, so we have to use
  // uint8_t here which is definitely 1 byte in size.
  // Also, C++ guarantees conversion between booleans and integers:
  // true <--> 1, false <--> 0.
  auto a = std::array<uint8_t, bit_stride>();
  auto t = std::array<uint64_t, byte_stride>();
  size_t dest_idx = offset;
  auto itr_finish = src.cbegin() + (src.size() / bit_stride) * bit_stride;

  for (auto itr = src.cbegin(); itr != itr_finish; itr += bit_stride) {
    std::copy(itr, itr + bit_stride, a.begin());
    std::memcpy(t.data(), a.data(), a.size());
    std::transform(t.cbegin(), t.cend(), dest.begin() + dest_idx,
                   [](auto e) { return (magic * e) >> 56; });
    dest_idx += byte_stride;
  }

  // For the remaining bits, process 8 bits at a time.
  for (auto itr = itr_finish; itr != src.cend(); itr += 8) {
    std::copy(itr, itr + 8, a.begin());
    std::memcpy(t.data(), a.data(), 8);
    dest[dest_idx++] = (magic * t[0]) >> 56;
  }

  return RTNType::Good;
}

inline auto SZ3::SPERR::unpack_booleans(std::vector<bool>& dest,
                            const void* src,
                            size_t src_len,
                            size_t src_offset) -> RTNType
{
  // For some reason, a strided unpack implementation similar to the striding
  // strategy used in `pack_booleans()` would result in slower runtime...
  // See Github issue #122.

  if (src == nullptr)
    return RTNType::Error;

  if (src_len < src_offset)
    return RTNType::WrongLength;

  const size_t num_of_bytes = src_len - src_offset;

  // `dest` needs to have enough space allocated, as the API specifies.
  assert(dest.size() >= num_of_bytes * 8);

  const uint8_t* src_ptr = static_cast<const uint8_t*>(src) + src_offset;
  const uint64_t magic = 0x8040201008040201;
  const uint64_t mask = 0x8080808080808080;

#ifndef OMP_UNPACK_BOOLEANS
  // Serial implementation
  //
  auto a = std::array<uint8_t, 8>();
  auto t = uint64_t{0};
  auto b8 = uint64_t{0};
  auto dest_itr = dest.begin();
  for (size_t byte_idx = 0; byte_idx < num_of_bytes; byte_idx++) {
    b8 = *(src_ptr + byte_idx);
    t = ((magic * b8) & mask) >> 7;
    std::memcpy(a.data(), &t, 8);
    std::copy(a.cbegin(), a.cend(), dest_itr);
    dest_itr += 8;
  }
#else
  // Parallel implementation
  //
  // Because in most implementations std::vector<bool> is stored as uint64_t
  // values, we parallel in strides of 64 bits, or 8 bytes.
  const size_t stride_size = 8;
  const size_t num_of_strides = num_of_bytes / stride_size;

#pragma omp parallel for
  for (size_t stride = 0; stride < num_of_strides; stride++) {
    uint8_t a[64]{0};
    for (size_t byte = 0; byte < stride_size; byte++) {
      const uint8_t* ptr = src_ptr + stride * stride_size + byte;
      const uint64_t t = ((magic * (*ptr)) & mask) >> 7;
      std::memcpy(a + byte * 8, &t, 8);
    }
    for (size_t i = 0; i < 64; i++)
      dest[stride * 64 + i] = a[i];
  }
  // This loop is at most 7 iterations, so not to worry about parallel anymore.
  for (size_t byte_idx = stride_size * num_of_strides; byte_idx < num_of_bytes; byte_idx++) {
    const uint8_t* ptr = src_ptr + byte_idx;
    const uint64_t t = ((magic * (*ptr)) & mask) >> 7;
    uint8_t a[8]{0};
    std::memcpy(a, &t, 8);
    for (size_t i = 0; i < 8; i++)
      dest[byte_idx * 8 + i] = a[i];
  }
#endif

  return RTNType::Good;
}

inline auto SZ3::SPERR::pack_8_booleans(std::array<bool, 8> src) -> uint8_t
{
  // It turns out that C++ doesn't specify bool to be one byte,
  // so to be safe we copy the content of src to array of uint8_t.
  auto bytes = std::array<uint8_t, 8>();
  std::copy(src.cbegin(), src.cend(), bytes.begin());
  const uint64_t magic = 0x8040201008040201;
  uint64_t t = 0;
  std::memcpy(&t, bytes.data(), 8);
  uint8_t dest = (magic * t) >> 56;
  return dest;
}

inline auto SZ3::SPERR::unpack_8_booleans(uint8_t src) -> std::array<bool, 8>
{
  const uint64_t magic = 0x8040201008040201;
  const uint64_t mask = 0x8080808080808080;
  uint64_t t = ((magic * src) & mask) >> 7;
  // It turns out that C++ doesn't specify bool to be one byte,
  // so to be safe we use an array of uint8_t.
  auto bytes = std::array<uint8_t, 8>();
  std::memcpy(bytes.data(), &t, 8);
  auto b8 = std::array<bool, 8>();
  std::copy(bytes.cbegin(), bytes.cend(), b8.begin());
  return b8;
}

inline auto SZ3::SPERR::read_n_bytes(std::string filename, size_t n_bytes) -> vec8_type
{
  auto buf = vec8_type();

  auto closer = [](std::FILE* f) { std::fclose(f); };  // bypass a compiler warning
  std::unique_ptr<std::FILE, decltype(closer)> fp(std::fopen(filename.data(), "rb"), closer);

  if (!fp)
    return buf;

  // POSIX systems require the size of a file to be specified, so
  // one can fseek to the end of the file.
  auto sk = std::fseek(fp.get(), 0, SEEK_END);
  assert(sk == 0);
  if (std::ftell(fp.get()) < n_bytes)
    return buf;

  std::rewind(fp.get());
  buf.resize(n_bytes);
  if (std::fread(buf.data(), 1, n_bytes, fp.get()) != n_bytes)
    buf.clear();

  return buf;
}

template <typename T>
inline auto SZ3::SPERR::read_whole_file(std::string filename) -> vec_type<T>
{
  auto buf = vec_type<T>();

  auto closer = [](std::FILE* f) { std::fclose(f); };  // bypass a compiler warning
  std::unique_ptr<std::FILE, decltype(closer)> fp(std::fopen(filename.data(), "rb"), closer);
  if (!fp)
    return buf;

  // POSIX systems require the size of a file to be specified, so
  // one can fseek to the end of the file.
  auto sk = std::fseek(fp.get(), 0, SEEK_END);
  assert(sk == 0);
  const size_t file_size = std::ftell(fp.get());
  if (file_size % sizeof(T) != 0)
    return buf;

  const size_t num_vals = file_size / sizeof(T);
  buf.resize(num_vals);
  std::rewind(fp.get());
  size_t nread = std::fread(buf.data(), sizeof(T), num_vals, fp.get());
  if (nread != num_vals)
    buf.clear();

  return buf;
}

inline auto SZ3::SPERR::write_n_bytes(std::string filename, size_t n_bytes, const void* buffer) -> RTNType
{
  auto closer = [](std::FILE* f) { std::fclose(f); };  // bypass a compiler warning
  std::unique_ptr<std::FILE, decltype(closer)> fp(std::fopen(filename.data(), "wb"), closer);
  if (!fp)
    return RTNType::IOError;

  if (std::fwrite(buffer, 1, n_bytes, fp.get()) != n_bytes)
    return RTNType::IOError;
  else
    return RTNType::Good;
}

inline auto SZ3::SPERR::read_sections(std::string filename,
                          const std::vector<size_t>& sections,
                          vec8_type& dst) -> RTNType
{
  // Calculate the farthest file location to be read.
  size_t far = 0;
  for (size_t i = 0; i < sections.size() / 2; i++)
    far = std::max(far, sections[i * 2] + sections[i * 2 + 1]);

  // Prepare to read the file.
  auto closer = [](std::FILE* f) { std::fclose(f); };  // bypass a compiler warning
  std::unique_ptr<std::FILE, decltype(closer)> fp(std::fopen(filename.data(), "rb"), closer);
  if (!fp)
    return RTNType::IOError;

  // Retrieve the file length in bytes.
  //    P.S. POSIX systems require the size of a file to be specified, so
  //    one can fseek to the end of the file.
  auto sk = std::fseek(fp.get(), 0, SEEK_END);
  assert(sk == 0);
  const size_t file_len = std::ftell(fp.get());
  if (file_len < far)
    return RTNType::WrongLength;

  // Calculate the resulting size of `dst`, and allocate enough memory.
  auto dst_pos = dst.size();  // keep track of the current position to write section data.
  auto total_len = dst.size();
  for (size_t i = 0; i < sections.size() / 2; i++)
    total_len += sections[i * 2 + 1];
  dst.resize(total_len);

  // Read in sections of the file!
  for (size_t i = 0; i < sections.size() / 2; i++) {
    sk = std::fseek(fp.get(), sections[i * 2], SEEK_SET);
    assert(sk == 0);
    auto nread = std::fread(dst.data() + dst_pos, 1, sections[i * 2 + 1], fp.get());
    assert(nread == sections[i * 2 + 1]);
    dst_pos += nread;
  }

  return RTNType::Good;
}

inline auto SZ3::SPERR::extract_sections(const void* buf,
                             size_t buf_len,
                             const std::vector<size_t>& sections,
                             vec8_type& dst) -> RTNType
{
  // Calculate the farthest file location to be read.
  size_t far = 0;
  for (size_t i = 0; i < sections.size() / 2; i++)
    far = std::max(far, sections[i * 2] + sections[i * 2 + 1]);
  if (buf_len < far)
    return RTNType::WrongLength;

  // Calculate the resulting size of `dst`, and allocate enough memory.
  auto total_len = dst.size();
  for (size_t i = 0; i < sections.size() / 2; i++)
    total_len += sections[i * 2 + 1];
  dst.reserve(total_len);

  // Extract sections of the buffer!
  for (size_t i = 0; i < sections.size() / 2; i++) {
    const auto* beg = static_cast<const uint8_t*>(buf) + sections[i * 2];
    const auto* end = beg + sections[i * 2 + 1];
    std::copy(beg, end, std::back_inserter(dst));
  }

  return RTNType::Good;
}

template <typename T>
inline auto SZ3::SPERR::calc_stats(const T* arr1, const T* arr2, size_t arr_len, size_t omp_nthreads)
    -> std::array<T, 5>
{
  const size_t stride_size = 8192;
  const size_t num_of_strides = arr_len / stride_size;
  const size_t remainder_size = arr_len - stride_size * num_of_strides;

  // Use the maximum possible threads if 0 is passed in.
#ifdef USE_OMP
  if (omp_nthreads == 0)
    omp_nthreads = omp_get_max_threads();
#endif

  auto rmse = T{0.0};
  auto linfty = T{0.0};
  auto psnr = T{0.0};
  auto arr1min = T{0.0};
  auto arr1max = T{0.0};

  //
  // Calculate min and max of arr1
  //
  const auto minmax = std::minmax_element(arr1, arr1 + arr_len);
  arr1min = *minmax.first;
  arr1max = *minmax.second;

  //
  // In rare cases, the two input arrays are identical.
  //
  auto is_equal = std::equal(arr1, arr1 + arr_len, arr2);
  if (is_equal) {
    psnr = std::numeric_limits<T>::infinity();
    return {rmse, linfty, psnr, arr1min, arr1max};
  }

  auto sum_vec = vec_type<T>(num_of_strides + 1);
  auto linfty_vec = vec_type<T>(num_of_strides + 1);

//
// Calculate diff summation and l-infty of each stride
//
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t stride_i = 0; stride_i < num_of_strides; stride_i++) {
    T maxerr = 0.0;
    auto buf = std::array<T, stride_size>();
    for (size_t i = 0; i < stride_size; i++) {
      const size_t idx = stride_i * stride_size + i;
      auto diff = std::abs(arr1[idx] - arr2[idx]);
      maxerr = std::max(maxerr, diff);
      buf[i] = diff * diff;
    }
    sum_vec[stride_i] = std::accumulate(buf.cbegin(), buf.cend(), T{0.0});
    linfty_vec[stride_i] = maxerr;
  }

  //
  // Calculate diff summation and l-infty of the remaining elements
  //
  T last_linfty = 0.0;
  auto last_buf = std::array<T, stride_size>{};  // must be enough for `remainder_size` elements.
  for (size_t i = 0; i < remainder_size; i++) {
    const size_t idx = stride_size * num_of_strides + i;
    auto diff = std::abs(arr1[idx] - arr2[idx]);
    last_linfty = std::max(last_linfty, diff);
    last_buf[i] = diff * diff;
  }
  sum_vec[num_of_strides] = T{0.0};
  sum_vec[num_of_strides] =
      std::accumulate(last_buf.cbegin(), last_buf.cbegin() + remainder_size, T{0.0});
  linfty_vec[num_of_strides] = last_linfty;

  //
  // Now calculate linfty
  //
  linfty = *(std::max_element(linfty_vec.cbegin(), linfty_vec.cend()));

  //
  // Now calculate rmse and psnr
  // Note: psnr is calculated in dB, and follows the equation described in:
  // http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/VELDHUIZEN/node18.html
  // Also refer to https://www.mathworks.com/help/vision/ref/psnr.html
  //
  const auto mse = std::accumulate(sum_vec.cbegin(), sum_vec.cend(), T{0.0}) / T(arr_len);
  rmse = std::sqrt(mse);
  const auto range_sq = (arr1max - arr1min) * (arr1max - arr1min);
  psnr = std::log10(range_sq / mse) * T{10.0};

  return {rmse, linfty, psnr, arr1min, arr1max};
}

template <typename T>
inline auto SZ3::SPERR::kahan_summation(const T* arr, size_t len) -> T
{
  T sum = 0.0, c = 0.0;
  T t, y;
  for (size_t i = 0; i < len; i++) {
    y = arr[i] - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }

  return sum;
}

inline auto SZ3::SPERR::chunk_volume(dims_type vol_dim,
                         dims_type chunk_dim) -> std::vector<std::array<size_t, 6>>
{
  // Step 1: figure out how many segments are there along each axis.
  auto n_segs = std::array<size_t, 3>();
  for (size_t i = 0; i < 3; i++) {
    n_segs[i] = vol_dim[i] / chunk_dim[i];
    // In case the last segment is shorter than `chunk_dim[i]`, test if it's
    // longer than half of `chunk_dim[i]`. If it is, it can qualify another segment.
    if ((vol_dim[i] % chunk_dim[i]) > (chunk_dim[i] / 2))
      n_segs[i]++;
    // In case the volume has one dimension that's too small, let's make it
    // at least one segment anyway.
    if (n_segs[i] == 0)
      n_segs[i] = 1;
  }

  // Step 2: calculate the starting indices of each segment along each axis.
  auto x_tics = std::vector<size_t>(n_segs[0] + 1);
  for (size_t i = 0; i < n_segs[0]; i++)
    x_tics[i] = i * chunk_dim[0];
  x_tics[n_segs[0]] = vol_dim[0];  // the last tic is the length in X

  auto y_tics = std::vector<size_t>(n_segs[1] + 1);
  for (size_t i = 0; i < n_segs[1]; i++)
    y_tics[i] = i * chunk_dim[1];
  y_tics[n_segs[1]] = vol_dim[1];  // the last tic is the length in Y

  auto z_tics = std::vector<size_t>(n_segs[2] + 1);
  for (size_t i = 0; i < n_segs[2]; i++)
    z_tics[i] = i * chunk_dim[2];
  z_tics[n_segs[2]] = vol_dim[2];  // the last tic is the length in Z

  // Step 3: fill in details of each chunk
  auto n_chunks = n_segs[0] * n_segs[1] * n_segs[2];
  auto chunks = std::vector<std::array<size_t, 6>>(n_chunks);
  size_t chunk_idx = 0;
  for (size_t z = 0; z < n_segs[2]; z++)
    for (size_t y = 0; y < n_segs[1]; y++)
      for (size_t x = 0; x < n_segs[0]; x++) {
        chunks[chunk_idx][0] = x_tics[x];                  // X start
        chunks[chunk_idx][1] = x_tics[x + 1] - x_tics[x];  // X length
        chunks[chunk_idx][2] = y_tics[y];                  // Y start
        chunks[chunk_idx][3] = y_tics[y + 1] - y_tics[y];  // Y length
        chunks[chunk_idx][4] = z_tics[z];                  // Z start
        chunks[chunk_idx][5] = z_tics[z + 1] - z_tics[z];  // Z length
        chunk_idx++;
      }

  return chunks;
}

template <typename T>
inline auto SZ3::SPERR::calc_mean_var(const T* arr, size_t len, size_t omp_nthreads) -> std::array<T, 2>
{
  if (len == 0) {
    static_assert(std::is_floating_point_v<T>);
    if constexpr (std::is_same_v<T, float>)
      return {std::nanf("1"), std::nanf("2")};
    else
      return {std::nan("1"), std::nan("2")};
  }

#ifdef USE_OMP
  if (omp_nthreads == 0)
    omp_nthreads = omp_get_max_threads();
#endif

  const size_t stride_size = 16'384;
  const size_t num_strides = len / stride_size;
  auto tmp_buf = vec_type<T>(num_strides + 1);

  // First, calculate the mean of this array.
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t i = 0; i < num_strides; i++) {
    const T* beg = arr + i * stride_size;
    tmp_buf[i] = std::accumulate(beg, beg + stride_size, T{0.0});
  }
  tmp_buf[num_strides] = 0.0;
  tmp_buf[num_strides] = std::accumulate(arr + num_strides * stride_size, arr + len, T{0.0});
  const auto elem_sum = std::accumulate(tmp_buf.cbegin(), tmp_buf.cend(), T{0.0});
  const auto mean = elem_sum / static_cast<T>(len);

  // Second, calculate the variance.
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t i = 0; i < num_strides; i++) {
    const T* beg = arr + i * stride_size;
    tmp_buf[i] = std::accumulate(beg, beg + stride_size, T{0.0}, [mean](auto init, auto v) {
      return init + (v - mean) * (v - mean);
    });
  }
  tmp_buf[num_strides] = 0.0;
  tmp_buf[num_strides] =
      std::accumulate(arr + num_strides * stride_size, arr + len, T{0.0},
                      [mean](auto init, auto v) { return init + (v - mean) * (v - mean); });
  const auto diff_sum = std::accumulate(tmp_buf.cbegin(), tmp_buf.cend(), T{0.0});
  const auto var = diff_sum / static_cast<T>(len);

  return {mean, var};
}

#ifdef __AVX2__
template <typename T>
inline auto SZ3::SPERR::any_ge_pow2(const T* buf, size_t len, T thld) -> bool
{
  assert((thld > 0) && (thld & (thld - 1)) == 0);

  const size_t simd_width = 32 / sizeof(T);
  T mask_val = ~(thld - 1);
  __m256i mask_vec;
  if constexpr (sizeof(T) == 8)
    mask_vec = _mm256_set1_epi64x(mask_val);
  else if constexpr (sizeof(T) == 4)
    mask_vec = _mm256_set1_epi32(mask_val);
  else if constexpr (sizeof(T) == 2)
    mask_vec = _mm256_set1_epi16(mask_val);
  else
    mask_vec = _mm256_set1_epi8(mask_val);

  size_t i = 0;
  for (; i + simd_width <= len; i += simd_width) {
    auto data_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(buf + i));
    if (!_mm256_testz_si256(data_vec, mask_vec))
      return true;
  }
  return std::any_of(buf + i, buf + len, [thld](auto v) { return v >= thld; });
}
#endif
// ---- END sperr_helper.cpp ----

// ---- BEGIN Bitstream.cpp ----
#include "Bitstream.h"

#include <cassert>
#include <cstring>
#include <iterator>  // std::distance()

#ifdef __SSE2__
#include <immintrin.h>
#endif

// Constructor
SZ3::SPERR::Bitstream::Bitstream(size_t nbits)
{
  m_itr = m_buf.begin();
  this->reserve(nbits);
}

// Functions for both read and write
inline void SZ3::SPERR::Bitstream::rewind()
{
  m_itr = m_buf.begin();
  m_buffer = 0;
  m_bits = 0;
}

inline auto SZ3::SPERR::Bitstream::capacity() const -> size_t
{
  return m_buf.size() * 64;
}

inline void SZ3::SPERR::Bitstream::reserve(size_t nbits)
{
  if (nbits > m_buf.size() * 64) {
#ifdef __SSE2__
    _mm_sfence();
#endif
    auto num_longs = (nbits + 63) / 64;
    const auto dist = std::distance(m_buf.begin(), m_itr);
    m_buf.resize(num_longs);         // trigger a memroy allocation.
    m_buf.resize(m_buf.capacity());  // be able to make use of all available capacity.
    m_itr = m_buf.begin() + dist;
  }
}

inline void SZ3::SPERR::Bitstream::reset()
{
  std::fill(m_buf.begin(), m_buf.end(), 0);
}

// Functions for read
inline auto SZ3::SPERR::Bitstream::rtell() const -> size_t
{
  // Stupid C++ insists that `m_buf.begin()` gives me a const iterator...
  std::vector<uint64_t>::const_iterator itr2 = m_itr;  // NOLINT
  return std::distance(m_buf.begin(), itr2) * 64 - m_bits;
}

inline void SZ3::SPERR::Bitstream::rseek(size_t offset)
{
  size_t div = offset >> 6;
  size_t rem = offset & 63;
  m_itr = m_buf.begin() + div;
  if (rem) {
    m_buffer = *m_itr >> rem;
    ++m_itr;
    m_bits = 64 - rem;
  }
  else {
    m_buffer = 0;
    m_bits = 0;
  }
}

inline auto SZ3::SPERR::Bitstream::rbit() -> bool
{
  if (m_bits == 0) {
    m_buffer = *m_itr;
    ++m_itr;
    m_bits = 64;
  }
  --m_bits;
  bool bit = m_buffer & uint64_t{1};
  m_buffer >>= 1;
  return bit;
}

// Functions for write
inline auto SZ3::SPERR::Bitstream::wtell() const -> size_t
{
  // Stupid C++ insists that `m_buf.begin()` gives me a const iterator...
  std::vector<uint64_t>::const_iterator itr2 = m_itr;  // NOLINT
  return std::distance(m_buf.begin(), itr2) * 64 + m_bits;
}

inline void SZ3::SPERR::Bitstream::wseek(size_t offset)
{
  size_t div = offset >> 6;
  size_t rem = offset & 63;
  m_itr = m_buf.begin() + div;
  if (rem) {
    m_buffer = *m_itr;
    m_buffer &= (uint64_t{1} << rem) - 1;
    m_bits = rem;
  }
  else {
    m_buffer = 0;
    m_bits = 0;
  }
}

inline void SZ3::SPERR::Bitstream::wbit(bool bit)
{
  m_buffer |= uint64_t{bit} << m_bits;

  if (++m_bits == 64) {
    if (m_itr == m_buf.end()) {  // allocate memory if necessary.
#ifdef __SSE2__
      _mm_sfence();
#endif
      auto dist = m_buf.size();
      m_buf.resize(std::max(size_t{1}, dist) * 2);
      m_itr = m_buf.begin() + dist;
    }
#ifdef __SSE2__
    auto dist = m_itr - m_buf.begin();
    long long int* ptr = reinterpret_cast<long long int*>(&m_buf[dist]);
    _mm_stream_si64(ptr, m_buffer);
#else
    *m_itr = m_buffer;
#endif
    ++m_itr;
    m_buffer = 0;
    m_bits = 0;
  }
}

inline void SZ3::SPERR::Bitstream::flush()
{
#ifdef __SSE2__
  _mm_sfence();
#endif
  if (m_bits) {  // only really flush when there are remaining bits.
    if (m_itr == m_buf.end()) {
      auto dist = m_buf.size();
      m_buf.resize(m_buf.size() + 1);
      m_itr = m_buf.begin() + dist;
    }
    *m_itr = m_buffer;
    ++m_itr;
    m_buffer = 0;
    m_bits = 0;
  }
}

// Functions to provide or parse a compact bitstream
inline void SZ3::SPERR::Bitstream::write_bitstream(void* p, size_t num_bits) const
{
  assert(num_bits <= m_buf.size() * 64);

  const auto num_longs = num_bits / 64;
  auto rem_bytes = num_bits / 8 - num_longs * sizeof(uint64_t);
  if (num_bits % 8 != 0)
    rem_bytes++;

  if (num_longs > 0)
    std::memcpy(p, m_buf.data(), num_longs * sizeof(uint64_t));

  if (rem_bytes > 0) {
    uint64_t value = m_buf[num_longs];
    auto* const p_byte = static_cast<std::byte*>(p);
    std::memcpy(p_byte + num_longs * sizeof(uint64_t), &value, rem_bytes);
  }
}

inline auto SZ3::SPERR::Bitstream::get_bitstream(size_t num_bits) const -> std::vector<std::byte>
{
  assert(num_bits <= m_buf.size() * 64);

  auto num_bytes = (num_bits + 7) / 8;
  auto tmp = std::vector<std::byte>(num_bytes);
  this->write_bitstream(tmp.data(), num_bits);

  return tmp;
}

inline void SZ3::SPERR::Bitstream::parse_bitstream(const void* p, size_t num_bits)
{
  this->reserve(num_bits);

  const auto num_longs = num_bits / 64;
  auto rem_bytes = num_bits / 8 - num_longs * sizeof(uint64_t);
  if (num_bits % 8 != 0)
    rem_bytes++;

  const auto* const p_byte = static_cast<const std::byte*>(p);

  if (num_longs > 0)
    std::memcpy(m_buf.data(), p_byte, num_longs * sizeof(uint64_t));

  if (rem_bytes > 0)
    std::memcpy(m_buf.data() + num_longs, p_byte + num_longs * sizeof(uint64_t), rem_bytes);

  this->rewind();
}
// ---- END Bitstream.cpp ----

// ---- BEGIN Bitmask.cpp ----
#include "Bitmask.h"

#include <algorithm>
#include <cassert>
#include <limits>

#if __cplusplus >= 202002L
#include <bit>
#endif

SZ3::SPERR::Bitmask::Bitmask(size_t nbits)
{
  auto num_longs = (nbits + 63) / 64;
  m_buf.assign(num_longs, 0);
  m_num_bits = nbits;
}

inline auto SZ3::SPERR::Bitmask::size() const -> size_t
{
  return m_num_bits;
}

inline void SZ3::SPERR::Bitmask::resize(size_t nbits)
{
  auto num_longs = (nbits + 63) / 64;
  m_buf.resize(num_longs, 0);
  m_num_bits = nbits;
}

inline auto SZ3::SPERR::Bitmask::rlong(size_t idx) const -> uint64_t
{
  return m_buf[idx >> 6];
}

inline auto SZ3::SPERR::Bitmask::rbit(size_t idx) const -> bool
{
  auto div = idx >> 6;  // idx / 64
  auto rem = idx & 63;  // idx % 64
  auto word = m_buf[div];
  word &= uint64_t{1} << rem;
  return word;
}

inline auto SZ3::SPERR::Bitmask::has_true(size_t start, size_t len) const -> bool
{
  auto long_idx = start >> 6;
  auto processed_bits = int64_t{0};
  auto word = m_buf[long_idx];

  // Collect the remaining bits from the start long.
  auto begin_idx = start & 63;
  auto nbits = std::min(size_t{64}, begin_idx + len);
  for (auto i = begin_idx; i < nbits; i++) {
    if (word & (uint64_t{1} << i))
      return true;
    processed_bits++;
  }

  // Examine the subsequent full longs.
  while (processed_bits + 64 <= len) {
    word = m_buf[++long_idx];
    if (word) {
      return true;
    }
    processed_bits += 64;
  }

  // Examine the remaining bits
  if (processed_bits < len) {
    nbits = len - processed_bits;
    assert(nbits < 64);
    word = m_buf[++long_idx];
    for (int64_t i = 0; i < nbits; i++) {
      if (word & (uint64_t{1} << i))
        return true;
    }
  }

  return false;
}

inline auto SZ3::SPERR::Bitmask::find_true(size_t start, size_t len) const -> int64_t
{
  auto long_idx = start >> 6;
  auto processed_bits = int64_t{0};
  auto word = m_buf[long_idx];

  // Collect the remaining bits from the start long.
  auto begin_idx = start & 63;
  auto nbits = std::min(size_t{64}, begin_idx + len);
  for (auto i = begin_idx; i < nbits; i++) {
    if (word & (uint64_t{1} << i))
      return processed_bits;
    processed_bits++;
  }

  // Examine the subsequent full longs.
  while (processed_bits + 64 <= len) {
    word = m_buf[++long_idx];
    if (word) {
#if __cplusplus >= 202002L
      int64_t i = std::countr_zero(word);
      return processed_bits + i;
#else
      for (int64_t i = 0; i < 64; i++)
        if (word & (uint64_t{1} << i))
          return processed_bits + i;
#endif
    }
    processed_bits += 64;
  }

  // Examine the remaining bits
  if (processed_bits < len) {
    nbits = len - processed_bits;
    assert(nbits < 64);
    word = m_buf[++long_idx];
    for (int64_t i = 0; i < nbits; i++) {
      if (word & (uint64_t{1} << i))
        return processed_bits + i;
    }
  }

  return -1;
}

inline auto SZ3::SPERR::Bitmask::count_true() const -> size_t
{
  size_t counter = 0;
  if (m_buf.empty())
    return counter;

  // Note that unused bits in the last long are not guaranteed to be all 0's.
  for (size_t i = 0; i < m_buf.size() - 1; i++) {
    auto val = m_buf[i];
#if __cplusplus >= 202002L
    counter += std::popcount(val);
#else
    if (val != 0) {
      for (size_t j = 0; j < 64; j++)
        counter += ((val >> j) & uint64_t{1});
    }
#endif
  }
  const auto val = m_buf.back();
  if (val != 0) {
    for (size_t j = 0; j < m_num_bits - (m_buf.size() - 1) * 64; j++)
      counter += ((val >> j) & uint64_t{1});
  }

  return counter;
}

inline void SZ3::SPERR::Bitmask::wlong(size_t idx, uint64_t value)
{
  m_buf[idx >> 6] = value;
}

inline void SZ3::SPERR::Bitmask::wbit(size_t idx, bool bit)
{
  const auto wstart = idx >> 6;
  auto word = m_buf[wstart];

  auto mask1 = uint64_t{1} << (idx & 63);
  word &= ~mask1;

  auto mask2 = uint64_t{bit} << (idx & 63);
  word |= mask2;

  m_buf[wstart] = word;
}

inline void SZ3::SPERR::Bitmask::wtrue(size_t idx)
{
  const auto wstart = idx >> 6;
  const auto mask = uint64_t{1} << (idx & 63);

  auto word = m_buf[wstart];
  word |= mask;
  m_buf[wstart] = word;
}

inline void SZ3::SPERR::Bitmask::wfalse(size_t idx)
{
  const auto wstart = idx >> 6;
  const auto mask = uint64_t{1} << (idx & 63);

  auto word = m_buf[wstart];
  word &= ~mask;
  m_buf[wstart] = word;
}

inline void SZ3::SPERR::Bitmask::reset()
{
  std::fill(m_buf.begin(), m_buf.end(), 0);
}

inline void SZ3::SPERR::Bitmask::reset_true()
{
  std::fill(m_buf.begin(), m_buf.end(), std::numeric_limits<uint64_t>::max());
}

inline auto SZ3::SPERR::Bitmask::view_buffer() const -> const std::vector<uint64_t>&
{
  return m_buf;
}

inline void SZ3::SPERR::Bitmask::use_bitstream(const void* p)
{
  const auto* pu64 = static_cast<const uint64_t*>(p);
  std::copy(pu64, pu64 + m_buf.size(), m_buf.begin());
}

#if __cplusplus >= 202002L && defined __cpp_lib_three_way_comparison
inline auto SZ3::SPERR::Bitmask::operator<=>(const Bitmask& rhs) const noexcept
{
  auto cmp = m_num_bits <=> rhs.m_num_bits;
  if (cmp != 0)
    return cmp;

  if (m_num_bits % 64 == 0)
    return m_buf <=> rhs.m_buf;
  else {
    // Compare each fully used long.
    for (size_t i = 0; i < m_buf.size() - 1; i++) {
      cmp = m_buf[i] <=> rhs.m_buf[i];
      if (cmp != 0)
        return cmp;
    }
    // Compare the last partially used long.
    auto mylast = m_buf.back();
    auto rhslast = rhs.m_buf.back();
    for (size_t i = m_num_bits % 64; i < 64; i++) {
      auto mask = uint64_t{1} << i;
      mylast &= ~mask;
      rhslast &= ~mask;
    }
    return mylast <=> rhslast;
  }
}
inline auto SZ3::SPERR::Bitmask::operator==(const Bitmask& rhs) const noexcept -> bool
{
  return (operator<=>(rhs) == 0);
}
#endif
// ---- END Bitmask.cpp ----

// ---- BEGIN Conditioner.cpp ----
#include <algorithm>  // std::all_of()
#include <cassert>
#include <cmath>    // std::sqrt()
#include <cstring>  // std::memcpy()
#include <numeric>  // std::accumulate()
#include <type_traits>

#include "Conditioner.h"

inline auto SZ3::SPERR::Conditioner::condition(vecd_type& buf, dims_type dims) -> condi_type
{
  // The order of performing condition operations:
  // 1. Test constant. If it's a constant field, return immediately.
  // 2. Subtract mean;

  assert(!buf.empty());
  auto meta = std::array<bool, 8>{true,    // subtract mean
                                  false,   // unused
                                  false,   // unused
                                  false,   // unused
                                  false,   // unused
                                  false,   // unused
                                  false,   // unused
                                  false};  // [7]: is this a constant field?

  // Operation 1
  //
  if (std::all_of(buf.cbegin(), buf.cend(), [v0 = buf[0]](auto v) { return v == v0; })) {
    meta[m_constant_field_idx] = true;
    const double val = buf[0];
    const uint64_t nval = buf.size();
    //
    // Assemble a header of the following info and ordering:
    // meta   nval  val
    //
    auto header = condi_type();
    header[0] = SZ3::SPERR::pack_8_booleans(meta);
    size_t pos = 1;
    std::memcpy(header.data() + pos, &nval, sizeof(nval));
    pos += sizeof(nval);
    std::memcpy(header.data() + pos, &val, sizeof(val));

    return header;
  }

  // Operation 2
  //
  m_adjust_strides(buf.size());
  const auto mean = m_calc_mean(buf);
  std::for_each(buf.begin(), buf.end(), [mean](auto& v) { v -= mean; });

  // Assemble a header of the following info order:
  // meta   mean  (empty)
  //
  auto header = condi_type();
  header[0] = SZ3::SPERR::pack_8_booleans(meta);
  size_t pos = 1;
  std::memcpy(header.data() + pos, &mean, sizeof(mean));
  pos += sizeof(mean);
  while (pos < header.size())
    header[pos++] = 0;

  return header;
}

inline auto SZ3::SPERR::Conditioner::inverse_condition(vecd_type& buf,
                                           dims_type dims,
                                           condi_type header) -> RTNType
{
  // unpack meta bit fields
  auto meta = SZ3::SPERR::unpack_8_booleans(header[0]);
  size_t pos = 1;

  // Operation 1: if this is a constant field?
  //
  if (meta[m_constant_field_idx]) {
    uint64_t nval = 0;
    double val = 0.0;
    std::memcpy(&nval, header.data() + pos, sizeof(nval));
    pos += sizeof(nval);
    std::memcpy(&val, header.data() + pos, sizeof(val));

    buf.resize(nval);
    std::fill(buf.begin(), buf.end(), val);

    return RTNType::Good;
  }

  // Operation 2: add back the mean
  //
  double mean = 0.0;
  std::memcpy(&mean, header.data() + pos, sizeof(mean));
  std::for_each(buf.begin(), buf.end(), [mean](auto& v) { v += mean; });

  return RTNType::Good;
}

inline auto SZ3::SPERR::Conditioner::is_constant(uint8_t byte) const -> bool
{
  auto b8 = SZ3::SPERR::unpack_8_booleans(byte);
  return b8[m_constant_field_idx];
}

inline void SZ3::SPERR::Conditioner::save_q(condi_type& header, double q) const
{
  // Save at position 9, the same as in `retrieve_q()`.
  std::memcpy(header.data() + 9, &q, sizeof(q));
}

inline auto SZ3::SPERR::Conditioner::retrieve_q(condi_type header) const -> double
{
  assert(!is_constant(header[0]));
  double q = 0.0;
  // Retrieve at position 9, the same as in `save_q()`.
  std::memcpy(&q, header.data() + 9, sizeof(q));
  return q;
}

inline auto SZ3::SPERR::Conditioner::m_calc_mean(const vecd_type& buf) -> double
{
  assert(buf.size() % m_num_strides == 0);

  m_stride_buf.resize(m_num_strides);
  const size_t stride_size = buf.size() / m_num_strides;

  for (size_t s = 0; s < m_num_strides; s++) {
    auto begin = buf.begin() + stride_size * s;
    auto end = begin + stride_size;
    m_stride_buf[s] = std::accumulate(begin, end, double{0.0}) / static_cast<double>(stride_size);
  }

  double sum = std::accumulate(m_stride_buf.begin(), m_stride_buf.end(), double{0.0});

  return (sum / static_cast<double>(m_stride_buf.size()));
}

inline void SZ3::SPERR::Conditioner::m_adjust_strides(size_t len)
{
  m_num_strides = m_default_num_strides;
  if (len % m_num_strides == 0)
    return;

  size_t num = 0;

  // First, try to increase till 2^15 = 32,768
  for (num = m_num_strides; num <= 32'768; num++) {
    if (len % num == 0)
      break;
  }

  if (len % num == 0) {
    m_num_strides = num;
    return;
  }

  // Second, try to decrease till 1, at which point it must work.
  for (num = m_num_strides; num > 0; num--) {
    if (len % num == 0)
      break;
  }

  m_num_strides = num;
}
// ---- END Conditioner.cpp ----

// ---- BEGIN CDF97.cpp ----
#include "CDF97.h"

#include <algorithm>
#include <cassert>
#include <numeric>  // std::accumulate()
#include <type_traits>

#ifdef __AVX2__
#include <immintrin.h>
#endif

// Destructor
SZ3::SPERR::CDF97::~CDF97()
{
  if (m_aligned_buf)
    SZ3::SPERR::aligned_free(m_aligned_buf);
}

template <typename T>
inline auto SZ3::SPERR::CDF97::copy_data(const T* data, size_t len, dims_type dims) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");
  if (len != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf.resize(len);
  std::copy(data, data + len, m_data_buf.begin());

  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * sizeof(double) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      SZ3::SPERR::aligned_free(m_aligned_buf);
    size_t alignment = 32;  // 256 bits
    size_t alloc_chunks = (max_col * 8 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<double*>(SZ3::SPERR::aligned_malloc(alignment, m_aligned_buf_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);

  return RTNType::Good;
}

inline auto SZ3::SPERR::CDF97::take_data(vecd_type&& buf, dims_type dims) -> RTNType
{
  if (buf.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf = std::move(buf);
  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * sizeof(double) > m_aligned_buf_bytes) {
    if (m_aligned_buf)
      SZ3::SPERR::aligned_free(m_aligned_buf);
    size_t alignment = 32;  // 256 bits
    size_t alloc_chunks = (max_col * 8 + 31) / alignment;
    m_aligned_buf_bytes = alignment * alloc_chunks;
    m_aligned_buf = static_cast<double*>(SZ3::SPERR::aligned_malloc(alignment, m_aligned_buf_bytes));
  }

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(max_slice);

  return RTNType::Good;
}

inline auto SZ3::SPERR::CDF97::view_data() const -> const vecd_type&
{
  return m_data_buf;
}

inline auto SZ3::SPERR::CDF97::release_data() -> vecd_type&&
{
  return std::move(m_data_buf);
}

inline auto SZ3::SPERR::CDF97::get_dims() const -> std::array<size_t, 3>
{
  return m_dims;
}

inline void SZ3::SPERR::CDF97::dwt1d()
{
  auto num_xforms = SZ3::SPERR::num_of_xforms(m_dims[0]);
  m_dwt1d(m_data_buf.data(), m_data_buf.size(), num_xforms);
}

inline void SZ3::SPERR::CDF97::idwt1d()
{
  auto num_xforms = SZ3::SPERR::num_of_xforms(m_dims[0]);
  m_idwt1d(m_data_buf.data(), m_data_buf.size(), num_xforms);
}

inline void SZ3::SPERR::CDF97::dwt2d()
{
  auto xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_dwt2d(m_data_buf.data(), {m_dims[0], m_dims[1]}, xy);
}

inline void SZ3::SPERR::CDF97::idwt2d()
{
  auto xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_idwt2d(m_data_buf.data(), {m_dims[0], m_dims[1]}, xy);
}

inline auto SZ3::SPERR::CDF97::idwt2d_multi_res() -> std::vector<vecd_type>
{
  const auto xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  auto ret = std::vector<vecd_type>();

  if (xy > 0) {
    ret.reserve(xy);
    for (size_t lev = xy; lev > 0; lev--) {
      auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(m_dims[1], lev);
      ret.emplace_back(m_sub_slice({x, y}));
      m_idwt2d_one_level(m_data_buf.data(), {x + xd, y + yd});
    }
  }

  return ret;
}

inline void SZ3::SPERR::CDF97::dwt3d()
{
  auto dyadic = SZ3::SPERR::can_use_dyadic(m_dims);
  if (dyadic)
    m_dwt3d_dyadic(*dyadic);
  else
    m_dwt3d_wavelet_packet();
}

inline void SZ3::SPERR::CDF97::idwt3d()
{
  auto dyadic = SZ3::SPERR::can_use_dyadic(m_dims);
  if (dyadic)
    m_idwt3d_dyadic(*dyadic);
  else
    m_idwt3d_wavelet_packet();
}

inline void SZ3::SPERR::CDF97::idwt3d_multi_res(std::vector<vecd_type>& h)
{
  auto dyadic = SZ3::SPERR::can_use_dyadic(m_dims);

  if (dyadic) {
    h.resize(*dyadic);
    for (size_t lev = *dyadic; lev > 0; lev--) {
      auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(m_dims[1], lev);
      auto [z, zd] = SZ3::SPERR::calc_approx_detail_len(m_dims[2], lev);
      auto& buf = h[*dyadic - lev];
      buf.resize(x * y * z);
      m_sub_volume({x, y, z}, buf.data());
      m_idwt3d_one_level({x + xd, y + yd, z + zd});
    }
  }
  else
    m_idwt3d_wavelet_packet();
}

inline void SZ3::SPERR::CDF97::m_dwt3d_wavelet_packet()
{
  /*
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First transform along the Z dimension
  //
  const auto num_xforms_z = SZ3::SPERR::num_of_xforms(m_dims[2]);

  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values of one XZ slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // DWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_dwt1d(m_slice_buf.data() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values of the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }

  // Second transform each plane
  //
  const auto num_xforms_xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));

  for (size_t z = 0; z < m_dims[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d(m_data_buf.data() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }
}

inline void SZ3::SPERR::CDF97::m_idwt3d_wavelet_packet()
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First, inverse transform each plane
  //
  auto num_xforms_xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  for (size_t i = 0; i < m_dims[2]; i++) {
    const size_t offset = plane_size_xy * i;
    m_idwt2d(m_data_buf.data() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }

  /*
   * Second, inverse transform along the Z dimension
   *
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  // Process one XZ slice at a time
  //
  const auto num_xforms_z = SZ3::SPERR::num_of_xforms(m_dims[2]);
  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values on one slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // IDWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_idwt1d(m_slice_buf.data() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values from the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }
}

inline void SZ3::SPERR::CDF97::m_dwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = 0; lev < num_xforms; lev++) {
    auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(m_dims[0], lev);
    auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(m_dims[1], lev);
    auto [z, zd] = SZ3::SPERR::calc_approx_detail_len(m_dims[2], lev);
    m_dwt3d_one_level({x, y, z});
  }
}

inline void SZ3::SPERR::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = SZ3::SPERR::calc_approx_detail_len(m_dims[2], lev - 1);
    m_idwt3d_one_level({x, y, z});
  }
}

//
// Private Methods
//
inline void SZ3::SPERR::CDF97::m_dwt1d(double* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    m_gather(array, array_len, m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, array_len);
    std::copy(m_aligned_buf, m_aligned_buf + array_len, array);
    array_len -= array_len / 2;
  }
}

inline void SZ3::SPERR::CDF97::m_idwt1d(double* array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(array_len, lev - 1);
    this->QccWAVCDF97SynthesisSymmetric(array, x);
    m_scatter(array, x, m_aligned_buf);
    std::copy(m_aligned_buf, m_aligned_buf + x, array);
  }
}

inline void SZ3::SPERR::CDF97::m_dwt2d(double* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(len_xy[0], lev);
    auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(len_xy[1], lev);
    m_dwt2d_one_level(plane, {x, y});
  }
}

inline void SZ3::SPERR::CDF97::m_idwt2d(double* plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = SZ3::SPERR::calc_approx_detail_len(len_xy[0], lev - 1);
    auto [y, yd] = SZ3::SPERR::calc_approx_detail_len(len_xy[1], lev - 1);
    m_idwt2d_one_level(plane, {x, y});
  }
}

inline void SZ3::SPERR::CDF97::m_dwt2d_one_level(double* plane, std::array<size_t, 2> len_xy)
{
  // First, perform DWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    m_gather(pos, len_xy[0], m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, len_xy[0]);
    std::copy(m_aligned_buf, m_aligned_buf + len_xy[0], pos);
  }

  // Second, perform DWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf[y] = plane[y * m_dims[0] + x];
    m_gather(m_slice_buf.data(), len_xy[1], m_aligned_buf);
    this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, len_xy[1]);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf[y];
  }
}

inline void SZ3::SPERR::CDF97::m_idwt2d_one_level(double* plane, std::array<size_t, 2> len_xy)
{
  // First, perform IDWT along Y for every column
  for (size_t x = 0; x < len_xy[0]; x++) {
    for (size_t y = 0; y < len_xy[1]; y++)
      m_slice_buf[y] = plane[y * m_dims[0] + x];
    this->QccWAVCDF97SynthesisSymmetric(m_slice_buf.data(), len_xy[1]);
    m_scatter(m_slice_buf.data(), len_xy[1], m_aligned_buf);
    for (size_t y = 0; y < len_xy[1]; y++)
      plane[y * m_dims[0] + x] = m_aligned_buf[y];
  }

  // Second, perform IDWT along X for every row
  for (size_t i = 0; i < len_xy[1]; i++) {
    auto* pos = plane + i * m_dims[0];
    this->QccWAVCDF97SynthesisSymmetric(pos, len_xy[0]);
    m_scatter(pos, len_xy[0], m_aligned_buf);
    std::copy(m_aligned_buf, m_aligned_buf + len_xy[0], pos);
  }
}

inline void SZ3::SPERR::CDF97::m_dwt3d_one_level(std::array<size_t, 3> len_xyz)
{
  // First, do one level of transform on all XY planes.
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto col_len = len_xyz[2];
  for (size_t z = 0; z < col_len; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d_one_level(m_data_buf.data() + offset, {len_xyz[0], len_xyz[1]});
  }

  // Second, do one level of transform on all Z columns.  Strategy:
  // 1) extract eight Z columns to buffer space `m_slice_buf`
  // 2) transform these eight columns
  // 3) put the Z columns back to their locations in the volume.
  //
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 8 doubles. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 7 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const size_t stride = std::min(size_t{8}, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf.data() + i * col_len;
        m_gather(itr, col_len, m_aligned_buf);
        this->QccWAVCDF97AnalysisSymmetric(m_aligned_buf, col_len);
        std::copy(m_aligned_buf, m_aligned_buf + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf[z + i * col_len];
      }
    }
  }
}

inline void SZ3::SPERR::CDF97::m_idwt3d_one_level(std::array<size_t, 3> len_xyz)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto col_len = len_xyz[2];

  // First, do one level of inverse transform on all Z columns.  Strategy:
  // 1) extract eight Z columns to buffer space `m_slice_buf`
  // 2) transform these eight columns
  // 3) put the Z columns back to their appropriate locations in the volume.
  //
  // Note: the reason to process eight columns at a time is that a cache line
  // is usually 64 bytes, or 8 doubles. That means when you pay the cost to retrieve
  // one value from the Z column, its neighboring 7 values are available for free!

  for (size_t y = 0; y < len_xyz[1]; y++) {
    for (size_t x = 0; x < len_xyz[0]; x += 8) {
      const size_t xy_offset = y * m_dims[0] + x;
      const size_t stride = std::min(size_t{8}, len_xyz[0] - x);

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_slice_buf[z + i * col_len] = m_data_buf[z * plane_size_xy + xy_offset + i];
      }

      for (size_t i = 0; i < stride; i++) {
        auto* itr = m_slice_buf.data() + i * col_len;
        this->QccWAVCDF97SynthesisSymmetric(itr, col_len);
        m_scatter(itr, col_len, m_aligned_buf);
        std::copy(m_aligned_buf, m_aligned_buf + col_len, itr);
      }

      for (size_t z = 0; z < col_len; z++) {
        for (size_t i = 0; i < stride; i++)
          m_data_buf[z * plane_size_xy + xy_offset + i] = m_slice_buf[z + i * col_len];
      }
    }
  }

  // Second, do one level of inverse transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level(m_data_buf.data() + offset, {len_xyz[0], len_xyz[1]});
  }
}

inline void SZ3::SPERR::CDF97::m_gather(const double* src, size_t len, double* dst) const
{
#ifdef __AVX2__
  const double* src_end = src + len;
  double* dst_evens = dst;
  double* dst_odds = dst + len - len / 2;

  // Process 8 elements at a time
  for (; src + 8 <= src_end; src += 8) {
    __m256d v0 = _mm256_loadu_pd(src);      // 0, 1, 2, 3
    __m256d v1 = _mm256_loadu_pd(src + 4);  // 4, 5, 6, 7

    __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
    __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

    __m256d result1 = _mm256_permute4x64_pd(evens, 0b11011000);  // 0, 2, 4, 6
    __m256d result2 = _mm256_permute4x64_pd(odds, 0b11011000);   // 1, 3, 5, 7

    _mm256_store_pd(dst_evens, result1);
    _mm256_storeu_pd(dst_odds, result2);

    dst_evens += 4;
    dst_odds += 4;
  }

  for (; src < src_end - 1; src += 2) {
    *(dst_evens++) = *src;
    *(dst_odds++) = *(src + 1);
  }

  if (src < src_end)
    *dst_evens = *src;
#else
  size_t low_count = len - len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dst = *(src + i * 2);
    ++dst;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dst = *(src + i * 2 + 1);
    ++dst;
  }
#endif
}

inline void SZ3::SPERR::CDF97::m_scatter(const double* begin, size_t len, double* dst) const
{
#ifdef __AVX2__
  const double* even_end = begin + len - len / 2;
  const double* odd_beg = even_end;
  const double* dst_end = dst + len;

  // Process 8 elements at a time
  for (; begin + 4 < even_end; begin += 4) {
    __m256d v0 = _mm256_loadu_pd(begin);    // 0, 1, 2, 3
    __m256d v1 = _mm256_loadu_pd(odd_beg);  // 4, 5, 6, 7

    __m256d evens = _mm256_unpacklo_pd(v0, v1);  // 0, 4, 2, 6
    __m256d odds = _mm256_unpackhi_pd(v0, v1);   // 1, 5, 3, 7

    __m256d result1 = _mm256_permute2f128_pd(evens, odds, 0x20);  // 0, 4, 1, 5
    __m256d result2 = _mm256_permute2f128_pd(evens, odds, 0x31);  // 2, 6, 3, 7

    _mm256_store_pd(dst, result1);
    _mm256_store_pd(dst + 4, result2);

    dst += 8;
    odd_beg += 4;
  }

  for (; dst < dst_end - 1; dst += 2) {
    *dst = *(begin++);
    *(dst + 1) = *(odd_beg++);
  }

  if (dst < dst_end)
    *dst = *begin;
#else
  size_t low_count = len - len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dst + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dst + i * 2 + 1) = *begin;
    ++begin;
  }
#endif
}

inline auto SZ3::SPERR::CDF97::m_sub_slice(std::array<size_t, 2> subdims) const -> vecd_type
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1]);

  auto ret = vecd_type(subdims[0] * subdims[1]);
  auto dst = ret.begin();
  for (size_t y = 0; y < subdims[1]; y++) {
    auto beg = m_data_buf.begin() + y * m_dims[0];
    std::copy(beg, beg + subdims[0], dst);
    dst += subdims[0];
  }

  return ret;
}

inline void SZ3::SPERR::CDF97::m_sub_volume(dims_type subdims, double* dst) const
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1] && subdims[2] <= m_dims[2]);

  const auto slice_len = m_dims[0] * m_dims[1];
  for (size_t z = 0; z < subdims[2]; z++) {
    for (size_t y = 0; y < subdims[1]; y++) {
      auto beg = m_data_buf.begin() + z * slice_len + y * m_dims[0];
      std::copy(beg, beg + subdims[0], dst);
      dst += subdims[0];
    }
  }
}

//
// Methods from QccPack
//
inline void SZ3::SPERR::CDF97::QccWAVCDF97AnalysisSymmetric(double* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  double* even = signal;
  double* odd = signal + even_len;

  // Process all the odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += ALPHA * (even[i] + even[i + 1]);
  odd[odd_len - 1] += ALPHA * (even[odd_len - 1] + even[even_len - 1]);

  // Process all the even elements
  even[0] += 2.0 * BETA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] += BETA * (odd[i - 1] + odd[i]);
  even[even_len - 1] += BETA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process all the odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] += GAMMA * (even[i] + even[i + 1]);
  odd[odd_len - 1] += GAMMA * (even[odd_len - 1] + even[even_len - 1]);

  // Process even elements
  even[0] = EPSILON * (even[0] + 2.0 * DELTA * odd[0]);
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = EPSILON * (even[i] + DELTA * (odd[i - 1] + odd[i]));
  even[even_len - 1] =
      EPSILON * (even[even_len - 1] + DELTA * (odd[even_len - 2] + odd[odd_len - 1]));

  // Process odd elements
  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= -INV_EPSILON;
}

inline void SZ3::SPERR::CDF97::QccWAVCDF97SynthesisSymmetric(double* signal, size_t len)
{
  size_t even_len = len - len / 2;
  size_t odd_len = len / 2;
  double* even = signal;
  double* odd = signal + even_len;

  // Process odd elements
  for (size_t i = 0; i < odd_len; i++)
    odd[i] *= (-EPSILON);

  // Process even elements
  even[0] = even[0] * INV_EPSILON - 2.0 * DELTA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] = even[i] * INV_EPSILON - DELTA * (odd[i - 1] + odd[i]);
  even[even_len - 1] =
      even[even_len - 1] * INV_EPSILON - DELTA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= GAMMA * (even[i] + even[i + 1]);
  odd[odd_len - 1] -= GAMMA * (even[odd_len - 1] + even[even_len - 1]);

  // Process even elements
  even[0] -= 2.0 * BETA * odd[0];
  for (size_t i = 1; i < even_len - 1; i++)
    even[i] -= BETA * (odd[i - 1] + odd[i]);
  even[even_len - 1] -= BETA * (odd[even_len - 2] + odd[odd_len - 1]);

  // Process odd elements
  for (size_t i = 0; i < odd_len - 1; i++)
    odd[i] -= ALPHA * (even[i] + even[i + 1]);
  odd[odd_len - 1] -= ALPHA * (even[odd_len - 1] + even[even_len - 1]);
}
// ---- END CDF97.cpp ----

// ---- BEGIN SPECK_INT.cpp ----
#include "SPECK_INT.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

#if __cplusplus >= 202002L
#include <bit>
#endif

//
// Free-standing helper function
//
inline auto SZ3::SPERR::speck_int_get_num_bitplanes(const void* buf) -> uint8_t
{
  // Given the header definition, directly retrieve the value stored in the first byte.
  const auto* const ptr = static_cast<const uint8_t*>(buf);
  return ptr[0];
}

template <typename T>
SZ3::SPERR::SPECK_INT<T>::SPECK_INT()
{
  static_assert(std::is_integral_v<T>);
  static_assert(std::is_unsigned_v<T>);
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::integer_len() const -> size_t
{
  if constexpr (std::is_same_v<uint64_t, T>)
    return sizeof(uint64_t);
  else if constexpr (std::is_same_v<uint32_t, T>)
    return sizeof(uint32_t);
  else if constexpr (std::is_same_v<uint16_t, T>)
    return sizeof(uint16_t);
  else
    return sizeof(uint8_t);
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::set_dims(dims_type dims)
{
  m_dims = dims;
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::set_budget(size_t bud)
{
  if (bud == 0)
    m_budget = std::numeric_limits<size_t>::max();
  else {
    while (bud % 8 != 0)
      bud++;
    m_budget = bud;
  }
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::get_speck_num_bits(const void* buf) const -> uint64_t
{
  // Given the header definition, directly retrieve the value stored in bytes 1--9.
  const auto* const ptr = static_cast<const uint8_t*>(buf);
  uint64_t num_bits = 0;
  std::memcpy(&num_bits, ptr + 1, sizeof(num_bits));
  return num_bits;
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::get_stream_full_len(const void* buf) const -> uint64_t
{
  auto num_bits = get_speck_num_bits(buf);
  while (num_bits % 8 != 0)
    ++num_bits;
  return (header_size + num_bits / 8);
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::use_bitstream(const void* p, size_t len)
{
  // Header definition: 9 bytes in total:
  // num_bitplanes (uint8_t), num_useful_bits (uint64_t)

  // Step 1: extract num_bitplanes and num_useful_bits
  assert(len >= header_size);
  const auto* const p8 = static_cast<const uint8_t*>(p);
  std::memcpy(&m_num_bitplanes, p8, sizeof(m_num_bitplanes));
  std::memcpy(&m_total_bits, p8 + sizeof(m_num_bitplanes), sizeof(m_total_bits));

  // Step 2: unpack bits.
  //    Note that the bitstream passed in might not be of its original length as a result of
  //    progressive access. In that case, we parse available bits, and pad 0's to make the
  //    bitstream still have `m_total_bits`.
  m_avail_bits = (len - header_size) * 8;
  if (m_avail_bits < m_total_bits) {
    m_bit_buffer.reserve(m_total_bits);
    m_bit_buffer.reset();  // Set buffer to contain all 0's.
    m_bit_buffer.parse_bitstream(p8 + header_size, m_avail_bits);
  }
  else {
    assert(m_avail_bits - m_total_bits < 64);
    m_avail_bits = m_total_bits;
    m_bit_buffer.parse_bitstream(p8 + header_size, m_total_bits);
  }

  // After parsing an incoming bitstream, m_avail_bits <= m_total_bits.
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::encode()
{
  m_initialize_lists();
  const auto coeff_len = m_dims[0] * m_dims[1] * m_dims[2];
  m_bit_buffer.reserve(coeff_len);  // A good starting point
  m_bit_buffer.rewind();
  m_total_bits = 0;

  // Mark every coefficient as insignificant
  m_LSP_mask.resize(coeff_len);
  m_LSP_mask.reset();
  m_LSP_new.clear();
  m_LSP_new.reserve(coeff_len / 16);
  m_LIP_mask.resize(coeff_len);
  m_LIP_mask.reset();

  // Treat it as a special case when all coeffs (m_coeff_buf) are zero.
  //    In such a case, we mark `m_num_bitplanes` as zero.
  //    Of course, `m_total_bits` is also zero.
  if (std::all_of(m_coeff_buf.cbegin(), m_coeff_buf.cend(), [](auto v) { return v == 0; })) {
    m_num_bitplanes = 0;
    return;
  }

  // Decide the starting threshold.
  const auto max_coeff = *std::max_element(m_coeff_buf.cbegin(), m_coeff_buf.cend());
  m_num_bitplanes = 1;
  m_threshold = 1;
  // !! Careful loop condition so no integer overflow !!
  while (max_coeff - m_threshold >= m_threshold) {
    m_threshold *= uint_type{2};
    m_num_bitplanes++;
  }

  // Marching over bitplanes.
  for (uint8_t bitplane = 0; bitplane < m_num_bitplanes; bitplane++) {
    m_sorting_pass();
    if (m_bit_buffer.wtell() >= m_budget)  // Happens only when fixed-rate compression.
      break;

    m_refinement_pass_encode();
    if (m_bit_buffer.wtell() >= m_budget)  // Happens only when fixed-rate compression.
      break;

    m_threshold /= uint_type{2};
    m_clean_LIS();
  }

  // Record the total number of bits produced, and flush the stream.
  m_total_bits = m_bit_buffer.wtell();
  m_bit_buffer.flush();
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::decode()
{
  m_initialize_lists();
  m_bit_buffer.rewind();

  // initialize coefficients to be zero, and sign array to be all positive
  const auto coeff_len = m_dims[0] * m_dims[1] * m_dims[2];
  m_coeff_buf.assign(coeff_len, uint_type{0});
  m_sign_array.resize(coeff_len);
  m_sign_array.reset_true();

  // Mark every coefficient as insignificant.
  m_LSP_mask.resize(coeff_len);
  m_LSP_mask.reset();
  m_LSP_new.clear();
  m_LSP_new.reserve(coeff_len / 16);
  m_LIP_mask.resize(coeff_len);
  m_LIP_mask.reset();

  // Handle the special case of all coeffs (m_coeff_buf) are zero by return now!
  // This case is indicated by both `m_num_bitplanes` and `m_total_bits` equal zero.
  if (m_num_bitplanes == 0) {
    assert(m_total_bits == 0);
    return;
  }

  // Restore the biggest `m_threshold`.
  m_threshold = 1;
  for (uint8_t i = 1; i < m_num_bitplanes; i++)
    m_threshold *= uint_type{2};

  // Marching over bitplanes.
  for (uint8_t bitplane = 0; bitplane < m_num_bitplanes; bitplane++) {
    m_sorting_pass();
    if (m_bit_buffer.rtell() >= m_avail_bits)  // Happens when a partial bitstream is available,
      break;                                   // because of progressive decoding or fixed-rate.

    m_refinement_pass_decode();
    if (m_bit_buffer.rtell() >= m_avail_bits)  // Happens when a partial bitstream is available,
      break;                                   // because of progressive decoding or fixed-rate.

    m_threshold /= uint_type{2};
    m_clean_LIS();
  }

  // The majority of newly identified significant points are initialized by the refinement pass.
  //    However, if the loop breaks after executing the sorting pass, then it leaves newly
  //    identified significant points from this iteration not initialized. We detect this case and
  //    initialize them here. The initialization strategy is the same as in the refinement pass.
  //
  if (!m_LSP_new.empty()) {
    const auto init_val = m_threshold + m_threshold - m_threshold / uint_type{2} - uint_type{1};
    for (auto idx : m_LSP_new)
      m_coeff_buf[idx] = init_val;
  }

  if (m_avail_bits == m_total_bits)
    assert(m_bit_buffer.rtell() == m_total_bits);
  else {
    assert(m_bit_buffer.rtell() >= m_avail_bits);
    assert(m_bit_buffer.rtell() <= m_total_bits);
  }
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::use_coeffs(vecui_type coeffs, Bitmask signs) -> RTNType
{
  if (coeffs.size() != signs.size())
    return RTNType::Error;
  m_coeff_buf = std::move(coeffs);
  m_sign_array = std::move(signs);
  return RTNType::Good;
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::release_coeffs() -> vecui_type&&
{
  return std::move(m_coeff_buf);
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::release_signs() -> Bitmask&&
{
  return std::move(m_sign_array);
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::view_coeffs() const -> const vecui_type&
{
  return m_coeff_buf;
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::view_signs() const -> const Bitmask&
{
  return m_sign_array;
}

template <typename T>
inline auto SZ3::SPERR::SPECK_INT<T>::encoded_bitstream_len() const -> size_t
{
  // Note that `m_total_bits` and `m_budget` can have 3 comparison outcomes:
  //  1. `m_total_bits < m_budget` no matter whether m_budget is the maximum size_t or not.
  //      In this case, we record all `m_total_bits` bits (simple).
  //  2. `m_total_bits > m_budget`: in this case, we record the value of `m_total_bits` in
  //      the header but only pack `m_budget` bits to the bitstream, creating an equivalence
  //      of truncating the first `m_budget` bits from a bitstream of length `m_total_bits`.
  //  3. `m_total_bits == m_budget`: this case is very unlikely, but if it happens, that's when
  //      `m_budget` happens to be exactly met after a sorting or refinement pass.
  //      In this case, we can also record all `m_total_bits` bits, same as outcome 1.
  //
  auto bits_to_pack = std::min(m_budget, size_t{m_total_bits});
  auto bit_in_byte = bits_to_pack / size_t{8};
  if (bits_to_pack % 8 != 0)
    ++bit_in_byte;
  return (header_size + bit_in_byte);
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::append_encoded_bitstream(vec8_type& buffer) const
{
  // Step 1: calculate size and allocate space for the encoded bitstream
  //
  // Header definition: 9 bytes in total:
  // num_bitplanes (uint8_t), num_useful_bits (uint64_t)
  //
  const auto app_size = this->encoded_bitstream_len();
  const auto orig_size = buffer.size();
  buffer.resize(orig_size + app_size);
  auto* const ptr = buffer.data() + orig_size;

  // Step 2: fill header
  size_t pos = 0;
  std::memcpy(ptr + pos, &m_num_bitplanes, sizeof(m_num_bitplanes));
  pos += sizeof(m_num_bitplanes);
  std::memcpy(ptr + pos, &m_total_bits, sizeof(m_total_bits));
  pos += sizeof(m_total_bits);

  // Step 3: assemble the right amount of bits into bytes.
  // See discussion on the number of bits to pack in function `encoded_bitstream_len()`.
  auto bits_to_pack = std::min(m_budget, size_t{m_total_bits});
  m_bit_buffer.write_bitstream(ptr + header_size, bits_to_pack);
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::m_refinement_pass_encode()
{
  // First, process significant pixels previously found.
  //
  const auto tmp1 = std::array<uint_type, 2>{uint_type{0}, m_threshold};
  const auto bits_x64 = m_LSP_mask.size() - m_LSP_mask.size() % 64;

  for (size_t i = 0; i < bits_x64; i += 64) {  // Evaluate 64 bits at a time.
    auto value = m_LSP_mask.rlong(i);

#if __cplusplus >= 202002L
    while (value) {
      auto j = std::countr_zero(value);
      const bool o1 = m_coeff_buf[i + j] >= m_threshold;
      m_coeff_buf[i + j] -= tmp1[o1];
      m_bit_buffer.wbit(o1);
      value &= value - 1;
    }
#else
    if (value != 0) {
      for (size_t j = 0; j < 64; j++) {
        if ((value >> j) & uint64_t{1}) {
          const bool o1 = m_coeff_buf[i + j] >= m_threshold;
          m_coeff_buf[i + j] -= tmp1[o1];
          m_bit_buffer.wbit(o1);
        }
      }
    }
#endif
  }
  for (auto i = bits_x64; i < m_LSP_mask.size(); i++) {  // Evaluate the remaining bits.
    if (m_LSP_mask.rbit(i)) {
      const bool o1 = m_coeff_buf[i] >= m_threshold;
      m_coeff_buf[i] -= tmp1[o1];
      m_bit_buffer.wbit(o1);
    }
  }

  // Second, mark newly found significant pixels in `m_LSP_mask`.
  //
  for (auto idx : m_LSP_new)
    m_LSP_mask.wtrue(idx);
  m_LSP_new.clear();
}

template <typename T>
inline void SZ3::SPERR::SPECK_INT<T>::m_refinement_pass_decode()
{
  // First, process significant points previously found.
  //    All these nested conditions are a little annoying, but I don't have a better solution now.
  //    Here's a documentation of their purposes.
  // 1) The decoding scheme (reconstructing values at the middle of an interval) requires
  //    different treatment when `m_threshold` is 1 or not.
  // 2) We make use of the internal representation of `m_LSP_mask` and evaluate 64 bits at time.
  //    This requires evaluating any remaining bits not divisible by 64.
  // 3) During progressive or fixed-rate decoding, we need to evaluate if the bitstream is
  //    exhausted after every read. We test it no matter what decoding mode we're in though.
  // 4) goto is used again. Here's it's used to jump out of a nested loop, which is an endorsed
  //    usage of it: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Res-goto
  //
  auto read_pos = m_bit_buffer.rtell();  // Avoid repeated calls to rtell().
  const auto bits_x64 = m_LSP_mask.size() - m_LSP_mask.size() % 64;  // <-- Point 2
  if (m_threshold >= uint_type{2}) {                                 // <-- Point 1
    const auto half_t = m_threshold / uint_type{2};
    for (size_t i = 0; i < bits_x64; i += 64) {  // <-- Point 2
      auto value = m_LSP_mask.rlong(i);

#if __cplusplus >= 202002L
      while (value) {
        auto j = std::countr_zero(value);
        if (m_bit_buffer.rbit())
          m_coeff_buf[i + j] += half_t;
        else
          m_coeff_buf[i + j] -= half_t;
        if (++read_pos == m_avail_bits)              // <-- Point 3
          goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;  // <-- Point 4
        value &= value - 1;
      }
#else
      if (value != 0) {
        for (size_t j = 0; j < 64; j++) {
          if ((value >> j) & uint64_t{1}) {
            if (m_bit_buffer.rbit())
              m_coeff_buf[i + j] += half_t;
            else
              m_coeff_buf[i + j] -= half_t;
            if (++read_pos == m_avail_bits)              // <-- Point 3
              goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;  // <-- Point 4
          }
        }
      }
#endif
    }
    for (auto i = bits_x64; i < m_LSP_mask.size(); i++) {  // <-- Point 2
      if (m_LSP_mask.rbit(i)) {
        if (m_bit_buffer.rbit())
          m_coeff_buf[i] += half_t;
        else
          m_coeff_buf[i] -= half_t;
        if (++read_pos == m_avail_bits)              // <-- Point 3
          goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;  // <-- Point 4
      }
    }
  }  // Finish the case where `m_threshold >= 2`.
  else {  // Start the case where `m_threshold == 1`.
    for (size_t i = 0; i < bits_x64; i += 64) {
      auto value = m_LSP_mask.rlong(i);

#if __cplusplus >= 202002L
      while (value) {
        auto j = std::countr_zero(value);
        if (m_bit_buffer.rbit())
          ++(m_coeff_buf[i + j]);
        if (++read_pos == m_avail_bits)
          goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;
        value &= value - 1;
      }
#else
      for (size_t j = 0; j < 64; j++) {
        if ((value >> j) & uint64_t{1}) {
          if (m_bit_buffer.rbit())
            ++(m_coeff_buf[i + j]);
          if (++read_pos == m_avail_bits)
            goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;
        }
      }
#endif
    }
    for (auto i = bits_x64; i < m_LSP_mask.size(); i++) {
      if (m_LSP_mask.rbit(i)) {
        if (m_bit_buffer.rbit())
          ++(m_coeff_buf[i]);
        if (++read_pos == m_avail_bits)
          goto INITIALIZE_NEWLY_FOUND_POINTS_LABEL;
      }
    }
  }
  assert(m_bit_buffer.rtell() <= m_avail_bits);

  // Second, initialize newly found significant points. Here I aim to initialize the reconstructed
  //    value at the middle of the interval specified by `m_threshold`. Note that given the integer
  //    nature of these coefficients, there are actually two values equally "in the middle."
  //    For example, with `m_threshold == 4`, the interval is [4, 8), and both 5 and 6 are "in the
  //    middle." I choose the smaller one (5 in this example) here. My experiments show that
  //    choosing the smaller value rather than the bigger one does not hurt, and sometimes bring a
  //    little extra PSNR gain (<0.5). Also note that the formula calculating `init_val`
  //    makes sure that when `m_threshold == 1`, significant points are initialized as 1.
  //
INITIALIZE_NEWLY_FOUND_POINTS_LABEL:
  const auto init_val = m_threshold + m_threshold - m_threshold / uint_type{2} - uint_type{1};
  for (auto idx : m_LSP_new)
    m_coeff_buf[idx] = init_val;
  for (auto idx : m_LSP_new)
    m_LSP_mask.wtrue(idx);
  m_LSP_new.clear();
}

// ---- END SPECK_INT.cpp ----

// ---- BEGIN SPECK3D_INT.cpp ----
#include "SPECK3D_INT.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

#if __cplusplus >= 202002L
#include <bit>
#endif

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT<T>::m_clean_LIS()
{
  for (auto& list : m_LIS) {
    auto it =
        std::remove_if(list.begin(), list.end(), [](const auto& s) { return s.num_elem() == 0; });
    list.erase(it, list.end());
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT<T>::m_initialize_lists()
{
  std::array<size_t, 3> num_of_parts;  // how many times each dimension could be partitioned?
  num_of_parts[0] = SZ3::SPERR::num_of_partitions(m_dims[0]);
  num_of_parts[1] = SZ3::SPERR::num_of_partitions(m_dims[1]);
  num_of_parts[2] = SZ3::SPERR::num_of_partitions(m_dims[2]);
  size_t num_of_sizes = std::accumulate(num_of_parts.cbegin(), num_of_parts.cend(), 1ul);

  // Initialize LIS
  if (m_LIS.size() < num_of_sizes)
    m_LIS.resize(num_of_sizes);
  std::for_each(m_LIS.begin(), m_LIS.end(), [](auto& list) { list.clear(); });

  // Starting from a set representing the whole volume, identify the smaller
  //    subsets and put them in the LIS accordingly.
  //    Note that it truncates 64-bit ints to 16-bit ints here, but should be OK.
  Set3D big;
  big.length_x = static_cast<uint16_t>(m_dims[0]);
  big.length_y = static_cast<uint16_t>(m_dims[1]);
  big.length_z = static_cast<uint16_t>(m_dims[2]);

  auto curr_lev = uint16_t{0};

  const auto dyadic = SZ3::SPERR::can_use_dyadic(m_dims);
  if (dyadic) {
    for (size_t i = 0; i < *dyadic; i++) {
      auto [subsets, next_lev] = m_partition_S_XYZ(big, curr_lev);
      big = subsets[0];
      for (auto it = std::next(subsets.cbegin()); it != subsets.cend(); ++it)
        m_LIS[next_lev].emplace_back(*it);
      curr_lev = next_lev;
    }
  }
  else {
    const auto num_xforms_xy = SZ3::SPERR::num_of_xforms(std::min(m_dims[0], m_dims[1]));
    const auto num_xforms_z = SZ3::SPERR::num_of_xforms(m_dims[2]);
    size_t xf = 0;
    while (xf < num_xforms_xy && xf < num_xforms_z) {
      auto [subsets, next_lev] = m_partition_S_XYZ(big, curr_lev);
      big = subsets[0];
      for (auto it = std::next(subsets.cbegin()); it != subsets.cend(); ++it)
        m_LIS[next_lev].emplace_back(*it);
      curr_lev = next_lev;
      xf++;
    }

    // One of these two conditions will happen.
    if (xf < num_xforms_xy) {
      while (xf < num_xforms_xy) {
        auto [subsets, next_lev] = m_partition_S_XY(big, curr_lev);
        big = subsets[0];
        for (auto it = std::next(subsets.cbegin()); it != subsets.cend(); ++it)
          m_LIS[next_lev].emplace_back(*it);
        curr_lev = next_lev;
        xf++;
      }
    }
    else if (xf < num_xforms_z) {
      while (xf < num_xforms_z) {
        auto [subsets, next_lev] = m_partition_S_Z(big, curr_lev);
        big = subsets[0];
        m_LIS[next_lev].emplace_back(subsets[1]);
        curr_lev = next_lev;
        xf++;
      }
    }
  }

  // Right now big is the set that's most likely to be significant, so insert
  // it at the front of it's corresponding vector. One-time expense.
  m_LIS[curr_lev].insert(m_LIS[curr_lev].begin(), big);

  // Encoder and decoder might have different additional tasks.
  m_additional_initialization();
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT<T>::m_sorting_pass()
{
  // Since we have a separate representation of LIP, let's process that list first!
  //
  const auto bits_x64 = m_LIP_mask.size() - m_LIP_mask.size() % 64;

  for (size_t i = 0; i < bits_x64; i += 64) {
    auto value = m_LIP_mask.rlong(i);

#if __cplusplus >= 202002L
    while (value) {
      auto j = std::countr_zero(value);
      m_process_P_lite(i + j);
      value &= value - 1;
    }
#else
    if (value != 0) {
      for (size_t j = 0; j < 64; j++) {
        if ((value >> j) & uint64_t{1})
          m_process_P_lite(i + j);
      }
    }
#endif
  }
  for (auto i = bits_x64; i < m_LIP_mask.size(); i++) {
    if (m_LIP_mask.rbit(i))
      m_process_P_lite(i);
  }

  // Then we process regular sets in LIS.
  //
  for (size_t tmp = 1; tmp <= m_LIS.size(); tmp++) {
    auto idx1 = m_LIS.size() - tmp;
    for (size_t idx2 = 0; idx2 < m_LIS[idx1].size(); idx2++) {
      size_t dummy = 0;
      m_process_S(idx1, idx2, dummy, true);
    }
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT<T>::m_code_S(size_t idx1, size_t idx2)
{
  auto set = m_LIS[idx1][idx2];

  if (set.length_x == 2 && set.length_y == 2 && set.length_z == 2) {  // tail ellison case
    size_t sig_counter = 0;
    bool need_decide = true;

    // Element (0, 0, 0)
    const auto id = set.start_z * m_dims[0] * m_dims[1] + set.start_y * m_dims[0] + set.start_x;
    auto mort = set.get_morton();
    m_LIP_mask.wtrue(id);
    m_process_P(id, mort, sig_counter, need_decide);

    // Element (1, 0, 0)
    auto id2 = id + 1;
    m_LIP_mask.wtrue(id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (0, 1, 0)
    id2 = id + m_dims[0];
    m_LIP_mask.wtrue(id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (1, 1, 0)
    m_LIP_mask.wtrue(++id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (0, 0, 1)
    id2 = id + m_dims[0] * m_dims[1];
    m_LIP_mask.wtrue(id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (1, 0, 1)
    m_LIP_mask.wtrue(++id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (0, 1, 1)
    id2 = id + m_dims[0] * (m_dims[1] + 1);
    m_LIP_mask.wtrue(id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);

    // Element (1, 1, 1)
    need_decide = sig_counter != 0;
    m_LIP_mask.wtrue(++id2);
    m_process_P(id2, ++mort, sig_counter, need_decide);
  }
  else {  // normal recursion case
          // Get its 8 subsets, and move the empty ones to the end.
    auto [subsets, next_lev] = m_partition_S_XYZ(set, uint16_t(idx1));
    const auto set_end =
        std::remove_if(subsets.begin(), subsets.end(), [](auto& s) { return s.num_elem() == 0; });

    // Counter for the number of discovered significant sets.
    //    If no significant subset is found yet, and we're already looking at the last subset,
    //    then we know that this last subset IS significant.
    size_t sig_counter = 0;
    for (auto it = subsets.begin(); it != set_end; ++it) {
      bool need_decide = (sig_counter != 0 || it + 1 != set_end);
      if (it->num_elem() == 1) {
        auto idx = it->start_z * m_dims[0] * m_dims[1] + it->start_y * m_dims[0] + it->start_x;
        m_LIP_mask.wtrue(idx);
        m_process_P(idx, it->get_morton(), sig_counter, need_decide);
      }
      else {
        m_LIS[next_lev].emplace_back(*it);
        const auto newidx2 = m_LIS[next_lev].size() - 1;
        m_process_S(next_lev, newidx2, sig_counter, need_decide);
      }
    }
  }
}

template <typename T>
inline auto SZ3::SPERR::SPECK3D_INT<T>::m_partition_S_XYZ(Set3D set, uint16_t lev) const
    -> std::tuple<std::array<Set3D, 8>, uint16_t>
{
  // Integer promotion rules (https://en.cppreference.com/w/c/language/conversion) say that types
  //    shorter than `int` are implicitly promoted to be `int` to perform calculations, so just
  //    keep them as `int` here because they'll involve in calculations later.
  //
  const auto split_x = std::array<int, 2>{set.length_x - set.length_x / 2, set.length_x / 2};
  const auto split_y = std::array<int, 2>{set.length_y - set.length_y / 2, set.length_y / 2};
  const auto split_z = std::array<int, 2>{set.length_z - set.length_z / 2, set.length_z / 2};

  const auto tmp = std::array<uint8_t, 2>{0, 1};
  lev += tmp[split_x[1] != 0];
  lev += tmp[split_y[1] != 0];
  lev += tmp[split_z[1] != 0];

  auto subsets = std::tuple<std::array<Set3D, 8>, uint16_t>();
  std::get<1>(subsets) = lev;
  auto morton_offset = set.get_morton();

  //
  // The actual figuring out where it starts/ends part...
  //
  // subset (0, 0, 0)
  auto& sub0 = std::get<0>(subsets)[0];
  sub0.set_morton(morton_offset);
  sub0.start_x = set.start_x;
  sub0.start_y = set.start_y;
  sub0.start_z = set.start_z;
  sub0.length_x = split_x[0];
  sub0.length_y = split_y[0];
  sub0.length_z = split_z[0];

  // subset (1, 0, 0)
  auto& sub1 = std::get<0>(subsets)[1];
  morton_offset += sub0.num_elem();
  sub1.set_morton(morton_offset);
  sub1.start_x = set.start_x + split_x[0];
  sub1.start_y = set.start_y;
  sub1.start_z = set.start_z;
  sub1.length_x = split_x[1];
  sub1.length_y = split_y[0];
  sub1.length_z = split_z[0];

  // subset (0, 1, 0)
  auto& sub2 = std::get<0>(subsets)[2];
  morton_offset += sub1.num_elem();
  sub2.set_morton(morton_offset);
  sub2.start_x = set.start_x;
  sub2.start_y = set.start_y + split_y[0];
  sub2.start_z = set.start_z;
  sub2.length_x = split_x[0];
  sub2.length_y = split_y[1];
  sub2.length_z = split_z[0];

  // subset (1, 1, 0)
  auto& sub3 = std::get<0>(subsets)[3];
  morton_offset += sub2.num_elem();
  sub3.set_morton(morton_offset);
  sub3.start_x = set.start_x + split_x[0];
  sub3.start_y = set.start_y + split_y[0];
  sub3.start_z = set.start_z;
  sub3.length_x = split_x[1];
  sub3.length_y = split_y[1];
  sub3.length_z = split_z[0];

  // subset (0, 0, 1)
  auto& sub4 = std::get<0>(subsets)[4];
  morton_offset += sub3.num_elem();
  sub4.set_morton(morton_offset);
  sub4.start_x = set.start_x;
  sub4.start_y = set.start_y;
  sub4.start_z = set.start_z + split_z[0];
  sub4.length_x = split_x[0];
  sub4.length_y = split_y[0];
  sub4.length_z = split_z[1];

  // subset (1, 0, 1)
  auto& sub5 = std::get<0>(subsets)[5];
  morton_offset += sub4.num_elem();
  sub5.set_morton(morton_offset);
  sub5.start_x = set.start_x + split_x[0];
  sub5.start_y = set.start_y;
  sub5.start_z = set.start_z + split_z[0];
  sub5.length_x = split_x[1];
  sub5.length_y = split_y[0];
  sub5.length_z = split_z[1];

  // subset (0, 1, 1)
  auto& sub6 = std::get<0>(subsets)[6];
  morton_offset += sub5.num_elem();
  sub6.set_morton(morton_offset);
  sub6.start_x = set.start_x;
  sub6.start_y = set.start_y + split_y[0];
  sub6.start_z = set.start_z + split_z[0];
  sub6.length_x = split_x[0];
  sub6.length_y = split_y[1];
  sub6.length_z = split_z[1];

  // subset (1, 1, 1)
  auto& sub7 = std::get<0>(subsets)[7];
  morton_offset += sub6.num_elem();
  sub7.set_morton(morton_offset);
  sub7.start_x = set.start_x + split_x[0];
  sub7.start_y = set.start_y + split_y[0];
  sub7.start_z = set.start_z + split_z[0];
  sub7.length_x = split_x[1];
  sub7.length_y = split_y[1];
  sub7.length_z = split_z[1];

  return subsets;
}

template <typename T>
inline auto SZ3::SPERR::SPECK3D_INT<T>::m_partition_S_XY(Set3D set, uint16_t lev) const
    -> std::tuple<std::array<Set3D, 4>, uint16_t>
{
  // This partition scheme is only used during initialization; no need to calculate morton offset.

  const auto split_x = std::array<int, 2>{set.length_x - set.length_x / 2, set.length_x / 2};
  const auto split_y = std::array<int, 2>{set.length_y - set.length_y / 2, set.length_y / 2};

  const auto tmp = std::array<uint8_t, 2>{0, 1};
  lev += tmp[split_x[1] != 0];
  lev += tmp[split_y[1] != 0];

  auto subsets = std::tuple<std::array<Set3D, 4>, uint16_t>();
  std::get<1>(subsets) = lev;
  const auto offsets = std::array<size_t, 3>{1, 2, 4};

  // The actual figuring out where it starts/ends part...
  //
  // subset (0, 0, 0)
  size_t sub_i = 0 * offsets[0] + 0 * offsets[1] + 0 * offsets[2];
  auto& sub0 = std::get<0>(subsets)[sub_i];
  sub0.start_x = set.start_x;
  sub0.start_y = set.start_y;
  sub0.start_z = set.start_z;
  sub0.length_x = split_x[0];
  sub0.length_y = split_y[0];
  sub0.length_z = set.length_z;

  // subset (1, 0, 0)
  sub_i = 1 * offsets[0] + 0 * offsets[1] + 0 * offsets[2];
  auto& sub1 = std::get<0>(subsets)[sub_i];
  sub1.start_x = set.start_x + split_x[0];
  sub1.start_y = set.start_y;
  sub1.start_z = set.start_z;
  sub1.length_x = split_x[1];
  sub1.length_y = split_y[0];
  sub1.length_z = set.length_z;

  // subset (0, 1, 0)
  sub_i = 0 * offsets[0] + 1 * offsets[1] + 0 * offsets[2];
  auto& sub2 = std::get<0>(subsets)[sub_i];
  sub2.start_x = set.start_x;
  sub2.start_y = set.start_y + split_y[0];
  sub2.start_z = set.start_z;
  sub2.length_x = split_x[0];
  sub2.length_y = split_y[1];
  sub2.length_z = set.length_z;

  // subset (1, 1, 0)
  sub_i = 1 * offsets[0] + 1 * offsets[1] + 0 * offsets[2];
  auto& sub3 = std::get<0>(subsets)[sub_i];
  sub3.start_x = set.start_x + split_x[0];
  sub3.start_y = set.start_y + split_y[0];
  sub3.start_z = set.start_z;
  sub3.length_x = split_x[1];
  sub3.length_y = split_y[1];
  sub3.length_z = set.length_z;

  return subsets;
}

template <typename T>
inline auto SZ3::SPERR::SPECK3D_INT<T>::m_partition_S_Z(Set3D set, uint16_t lev) const
    -> std::tuple<std::array<Set3D, 2>, uint16_t>
{
  // This partition scheme is only used during initialization; no need to calculate morton offset.

  const auto split_z = std::array<int, 2>{set.length_z - set.length_z / 2, set.length_z / 2};
  if (split_z[1] != 0)
    lev++;

  auto subsets = std::tuple<std::array<Set3D, 2>, uint16_t>();
  std::get<1>(subsets) = lev;

  //
  // The actual figuring out where it starts/ends part...
  //
  // subset (0, 0, 0)
  auto& sub0 = std::get<0>(subsets)[0];
  sub0.start_x = set.start_x;
  sub0.length_x = set.length_x;
  sub0.start_y = set.start_y;
  sub0.length_y = set.length_y;
  sub0.start_z = set.start_z;
  sub0.length_z = split_z[0];

  // subset (0, 0, 1)
  auto& sub1 = std::get<0>(subsets)[1];
  sub1.start_x = set.start_x;
  sub1.length_x = set.length_x;
  sub1.start_y = set.start_y;
  sub1.length_y = set.length_y;
  sub1.start_z = set.start_z + split_z[0];
  sub1.length_z = split_z[1];

  return subsets;
}

// ---- END SPECK3D_INT.cpp ----

// ---- BEGIN SPECK3D_INT_ENC.cpp ----
#include "SPECK3D_INT_ENC.h"

#include <algorithm>
#include <cassert>
#include <cstring>  // std::memcpy()
#include <numeric>

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_ENC<T>::m_deposit_set(Set3D set)
{
  switch (set.num_elem()) {
    case 0:
      return;
    case 1: {
      auto id = set.start_z * m_dims[0] * m_dims[1] + set.start_y * m_dims[0] + set.start_x;
      m_morton_buf[set.get_morton()] = m_coeff_buf[id];
      return;
    }
    case 2: {
      // We directly deposit the 2 elements in `set` instead of performing another partition.
      //
      // Deposit the 1st element.
      auto id = set.start_z * m_dims[0] * m_dims[1] + set.start_y * m_dims[0] + set.start_x;
      auto morton_id = set.get_morton();
      m_morton_buf[morton_id] = m_coeff_buf[id];

      // Deposit the 2nd element.
      if (set.length_x == 2)
        id++;
      else if (set.length_y == 2)
        id += m_dims[0];
      else
        id += m_dims[0] * m_dims[1];
      m_morton_buf[++morton_id] = m_coeff_buf[id];

      return;
    }
    case 4: {
      const auto id = set.start_z * m_dims[0] * m_dims[1] + set.start_y * m_dims[0] + set.start_x;
      auto morton_id = set.get_morton();

      if (set.length_x == 2 && set.length_y == 2) {
        // Element (0, 0, 0)
        m_morton_buf[morton_id] = m_coeff_buf[id];

        // Element (1, 0, 0)
        m_morton_buf[++morton_id] = m_coeff_buf[id + 1];

        // Element (0, 1, 0)
        auto id2 = id + m_dims[0];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (1, 1, 0)
        m_morton_buf[++morton_id] = m_coeff_buf[++id2];

        return;
      }
      else if (set.length_x == 2 && set.length_z == 2) {
        // Element (0, 0, 0)
        m_morton_buf[morton_id] = m_coeff_buf[id];

        // Element (1, 0, 0)
        m_morton_buf[++morton_id] = m_coeff_buf[id + 1];

        // Element (0, 0, 1)
        auto id2 = id + m_dims[0] * m_dims[1];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (1, 0, 1)
        m_morton_buf[++morton_id] = m_coeff_buf[++id2];

        return;
      }
      else if (set.length_y == 2 && set.length_z == 2) {
        // Element (0, 0, 0)
        m_morton_buf[morton_id] = m_coeff_buf[id];

        // Element (0, 1, 0)
        auto id2 = id + m_dims[0];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (0, 0, 1)
        id2 = id + m_dims[0] * m_dims[1];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (0, 1, 1)
        id2 += m_dims[0];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        return;
      }
      else
        break;  // Fall back to the recursive case.
    }
    case 8: {
      if (set.length_x == 2 && set.length_y == 2) {
        // Element (0, 0, 0)
        const auto id = set.start_z * m_dims[0] * m_dims[1] + set.start_y * m_dims[0] + set.start_x;
        auto morton_id = set.get_morton();
        m_morton_buf[morton_id] = m_coeff_buf[id];

        // Element (1, 0, 0)
        m_morton_buf[++morton_id] = m_coeff_buf[id + 1];

        // Element (0, 1, 0)
        auto id2 = id + m_dims[0];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (1, 1, 0)
        m_morton_buf[++morton_id] = m_coeff_buf[++id2];

        // Element (0, 0, 1)
        id2 = id + m_dims[0] * m_dims[1];
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (1, 0, 1)
        m_morton_buf[++morton_id] = m_coeff_buf[++id2];

        // Element (0, 1, 1)
        id2 = id + m_dims[0] * (m_dims[1] + 1);
        m_morton_buf[++morton_id] = m_coeff_buf[id2];

        // Element (1, 1, 1)
        m_morton_buf[++morton_id] = m_coeff_buf[++id2];

        return;
      }
      else
        break;  // Fall back to the recursive case.
    }
    default:
      break;  // Fall back to the recursive case.
  }

  // The recursive case.
  auto [subsets, lev] = m_partition_S_XYZ(set, 0);
  for (auto& sub : subsets)
    m_deposit_set(sub);
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_ENC<T>::m_additional_initialization()
{
  // For the encoder, this function re-organizes the coefficients in a morton order.
  //
  m_morton_buf.resize(m_dims[0] * m_dims[1] * m_dims[2]);

  // The same traversing order as in `SPECK3D_INT::m_sorting_pass()`
  size_t morton_offset = 0;
  for (size_t tmp = 1; tmp <= m_LIS.size(); tmp++) {
    auto idx1 = m_LIS.size() - tmp;
    for (size_t idx2 = 0; idx2 < m_LIS[idx1].size(); idx2++) {
      auto& set = m_LIS[idx1][idx2];
      set.set_morton(morton_offset);
      m_deposit_set(set);
      morton_offset += set.num_elem();
    }
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_ENC<T>::m_process_S(size_t idx1, size_t idx2, size_t& counter, bool output)
{
  auto& set = m_LIS[idx1][idx2];
  auto is_sig = true;

  // If need to output, it means the current set has unknown significance.
  if (output) {
    auto first = m_morton_buf.data() + set.get_morton();
#ifdef __AVX2__
    is_sig = SZ3::SPERR::any_ge_pow2(first, set.num_elem(), m_threshold);
#else
    is_sig = std::any_of(first, first + set.num_elem(),
                         [thld = m_threshold](auto v) { return v >= thld; });
#endif
    m_bit_buffer.wbit(is_sig);
  }

  if (is_sig) {
    counter++;
    m_code_S(idx1, idx2);
    set.make_empty();  // this current set is gonna be discarded.
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_ENC<T>::m_process_P(size_t idx, size_t morton, size_t& counter, bool output)
{
  bool is_sig = true;

  if (output) {
    assert(m_coeff_buf[idx] == m_morton_buf[morton]);
    is_sig = (m_morton_buf[morton] >= m_threshold);
    m_bit_buffer.wbit(is_sig);
  }

  if (is_sig) {
    counter++;  // Let's increment the counter first!
    assert(m_coeff_buf[idx] >= m_threshold);
    m_coeff_buf[idx] -= m_threshold;

    m_bit_buffer.wbit(m_sign_array.rbit(idx));
    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_ENC<T>::m_process_P_lite(size_t idx)
{
  auto is_sig = (m_coeff_buf[idx] >= m_threshold);
  m_bit_buffer.wbit(is_sig);

  if (is_sig) {
    assert(m_coeff_buf[idx] >= m_threshold);
    m_coeff_buf[idx] -= m_threshold;

    m_bit_buffer.wbit(m_sign_array.rbit(idx));
    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

// ---- END SPECK3D_INT_ENC.cpp ----

// ---- BEGIN SPECK3D_INT_DEC.cpp ----
#include "SPECK3D_INT_DEC.h"

#include <algorithm>
#include <cassert>
#include <cstring>  // std::memcpy()
#include <numeric>

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_DEC<T>::m_process_S(size_t idx1, size_t idx2, size_t& counter, bool read)
{
  auto& set = m_LIS[idx1][idx2];

  bool is_sig = true;
  if (read)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;
    m_code_S(idx1, idx2);
    set.make_empty();  // this current set is gonna be discarded.
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_DEC<T>::m_process_P(size_t idx, size_t no_use, size_t& counter, bool read)
{
  bool is_sig = true;
  if (read)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;  // Let's increment the counter first!
    m_sign_array.wbit(idx, m_bit_buffer.rbit());
    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK3D_INT_DEC<T>::m_process_P_lite(size_t idx)
{
  auto is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    m_sign_array.wbit(idx, m_bit_buffer.rbit());
    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

// ---- END SPECK3D_INT_DEC.cpp ----

// ---- BEGIN SPECK1D_INT.cpp ----
#include "SPECK1D_INT.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT<T>::m_clean_LIS()
{
  for (auto& list : m_LIS) {
    auto it =
        std::remove_if(list.begin(), list.end(), [](const auto& s) { return s.get_length() == 0; });
    list.erase(it, list.end());
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT<T>::m_initialize_lists()
{
  const auto total_len = m_dims[0];
  auto num_of_parts = SZ3::SPERR::num_of_partitions(total_len);
  auto num_of_lists = num_of_parts + 1;
  if (m_LIS.size() < num_of_lists)
    m_LIS.resize(num_of_lists);
  std::for_each(m_LIS.begin(), m_LIS.end(), [](auto& list) { list.clear(); });

  // Put in two sets, each representing a half of the long array.
  Set1D set;
  set.set_length(total_len);  // Set represents the whole 1D array.
  auto sets = m_partition_set(set);
  m_LIS[sets[0].get_level()].emplace_back(sets[0]);
  m_LIS[sets[1].get_level()].emplace_back(sets[1]);
}

template <typename T>
inline auto SZ3::SPERR::SPECK1D_INT<T>::m_partition_set(Set1D set) const -> std::array<Set1D, 2>
{
  const auto start = set.get_start();
  const auto length = set.get_length();
  const auto level = set.get_level();
  std::array<Set1D, 2> subsets;

  // Prepare the 1st set
  auto& set1 = subsets[0];
  set1.set_start(start);
  set1.set_length(length - length / 2);
  set1.set_level(level + 1);
  // Prepare the 2nd set
  auto& set2 = subsets[1];
  set2.set_start(start + length - length / 2);
  set2.set_length(length / 2);
  set2.set_level(level + 1);

  return subsets;
}

// ---- END SPECK1D_INT.cpp ----

// ---- BEGIN SPECK1D_INT_ENC.cpp ----
#include "SPECK1D_INT_ENC.h"

#include <algorithm>
#include <cassert>
#include <cstring>  // std::memcpy()
#include <numeric>

#if __cplusplus >= 202002L
#include <bit>
#endif

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_ENC<T>::m_sorting_pass()
{
  // Since we have a separate representation of LIP, let's process that list first!
  //
  const auto bits_x64 = m_LIP_mask.size() - m_LIP_mask.size() % 64;

  for (size_t i = 0; i < bits_x64; i += 64) {
    auto value = m_LIP_mask.rlong(i);

#if __cplusplus >= 202002L
    while (value) {
      size_t j = std::countr_zero(value);
      m_process_P(i + j, SigType::Dunno, j, true);
      value &= value - 1;
    }
#else
    if (value != 0) {
      for (size_t j = 0; j < 64; j++) {
        if ((value >> j) & uint64_t{1}) {
          size_t dummy = 0;
          m_process_P(i + j, SigType::Dunno, dummy, true);
        }
      }
    }
#endif
  }
  for (auto i = bits_x64; i < m_LIP_mask.size(); i++) {
    if (m_LIP_mask.rbit(i)) {
      size_t dummy = 0;
      m_process_P(i, SigType::Dunno, dummy, true);
    }
  }

  // Then we process regular sets in LIS.
  //
  for (size_t tmp = 1; tmp <= m_LIS.size(); tmp++) {
    // From the end of m_LIS to its front
    size_t idx1 = m_LIS.size() - tmp;
    for (size_t idx2 = 0; idx2 < m_LIS[idx1].size(); idx2++) {
      size_t dummy = 0;
      m_process_S(idx1, idx2, SigType::Dunno, dummy, true);
    }
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_ENC<T>::m_process_S(size_t idx1,
                                            size_t idx2,
                                            SigType sig,
                                            size_t& counter,
                                            bool output)
{
  auto& set = m_LIS[idx1][idx2];

  // Strategy to decide the significance of this set;
  // 1) If sig == dunno, then find the significance of this set. We do it in a
  //    way that at least one of its 2 subsets' significance become known as well.
  // 2) If sig is significant, then we directly proceed to `m_code_s()`, with its
  //    subsets' significance is unknown.
  // 3) if sig is insignificant, then this set is not processed.
  //
  auto subset_sigs = std::array<SigType, 2>{SigType::Dunno, SigType::Dunno};

  if (sig == SigType::Dunno) {
    auto set_sig = m_decide_significance(set);
    sig = set_sig ? SigType::Sig : SigType::Insig;
    if (set_sig) {
      if (*set_sig < set.get_length() - set.get_length() / 2)
        subset_sigs = {SigType::Sig, SigType::Dunno};
      else
        subset_sigs = {SigType::Insig, SigType::Sig};
    }
  }

  if (output)
    m_bit_buffer.wbit(sig == SigType::Sig);

  if (sig == SigType::Sig) {
    counter++;  // Let's increment the counter first!
    m_code_S(idx1, idx2, subset_sigs);
    set.set_length(0);  // this current set is gonna be discarded.
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_ENC<T>::m_process_P(size_t idx, SigType sig, size_t& counter, bool output)
{
  // Decide the significance of this pixel
  bool is_sig = false;
  if (sig == SigType::Dunno)
    is_sig = (m_coeff_buf[idx] >= m_threshold);
  else
    is_sig = (sig == SigType::Sig);

  if (output)
    m_bit_buffer.wbit(is_sig);

  if (is_sig) {
    counter++;  // Let's increment the counter first!
    m_bit_buffer.wbit(m_sign_array.rbit(idx));

    assert(m_coeff_buf[idx] >= m_threshold);
    m_coeff_buf[idx] -= m_threshold;
    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_ENC<T>::m_code_S(size_t idx1,
                                         size_t idx2,
                                         std::array<SigType, 2> subset_sigs)
{
  auto subsets = m_partition_set(m_LIS[idx1][idx2]);
  auto sig_counter = size_t{0};
  auto output = bool{true};

  // Process the 1st subset
  const auto& set0 = subsets[0];
  assert(set0.get_length() != 0);
  if (set0.get_length() == 1) {
    m_LIP_mask.wtrue(set0.get_start());
    m_process_P(set0.get_start(), subset_sigs[0], sig_counter, output);
  }
  else {
    const auto newidx1 = set0.get_level();
    m_LIS[newidx1].emplace_back(set0);
    const auto newidx2 = m_LIS[newidx1].size() - 1;
    m_process_S(newidx1, newidx2, subset_sigs[0], sig_counter, output);
  }

  // Process the 2nd subset
  if (sig_counter == 0) {
    output = false;
    subset_sigs[1] = SigType::Sig;
  }
  const auto& set1 = subsets[1];
  assert(set1.get_length() != 0);
  if (set1.get_length() == 1) {
    m_LIP_mask.wtrue(set1.get_start());
    m_process_P(set1.get_start(), subset_sigs[1], sig_counter, output);
  }
  else {
    const auto newidx1 = set1.get_level();
    m_LIS[newidx1].emplace_back(set1);
    const auto newidx2 = m_LIS[newidx1].size() - 1;
    m_process_S(newidx1, newidx2, subset_sigs[1], sig_counter, output);
  }
}

template <typename T>
inline auto SZ3::SPERR::SPECK1D_INT_ENC<T>::m_decide_significance(const Set1D& set) const
    -> std::optional<size_t>
{
  assert(set.get_length() != 0);

  const auto gtr = [thld = m_threshold](auto v) { return v >= thld; };
  auto first = m_coeff_buf.cbegin() + set.get_start();
  auto last = first + set.get_length();
  auto found = std::find_if(first, last, gtr);
  if (found != last)
    return static_cast<size_t>(std::distance(first, found));
  else
    return {};
}

// ---- END SPECK1D_INT_ENC.cpp ----

// ---- BEGIN SPECK1D_INT_DEC.cpp ----
#include "SPECK1D_INT_DEC.h"

#include <algorithm>
#include <cassert>
#include <cstring>  // std::memcpy()
#include <numeric>

#if __cplusplus >= 202002L
#include <bit>
#endif

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_DEC<T>::m_sorting_pass()
{
  // Since we have a separate representation of LIP, let's process that list first
  //
  const auto bits_x64 = m_LIP_mask.size() - m_LIP_mask.size() % 64;

  for (size_t i = 0; i < bits_x64; i += 64) {
    auto value = m_LIP_mask.rlong(i);

#if __cplusplus >= 202002L
    while (value) {
      size_t j = std::countr_zero(value);
      m_process_P(i + j, j, true);
      value &= value - 1;
    }
#else
    if (value != 0) {
      for (size_t j = 0; j < 64; j++) {
        if ((value >> j) & uint64_t{1}) {
          size_t dummy = 0;
          m_process_P(i + j, dummy, true);
        }
      }
    }
#endif
  }
  for (auto i = bits_x64; i < m_LIP_mask.size(); i++) {
    if (m_LIP_mask.rbit(i)) {
      size_t dummy = 0;
      m_process_P(i, dummy, true);
    }
  }

  // Then we process regular sets in LIS.
  //
  for (size_t tmp = 1; tmp <= m_LIS.size(); tmp++) {
    // From the end of m_LIS to its front
    size_t idx1 = m_LIS.size() - tmp;
    for (size_t idx2 = 0; idx2 < m_LIS[idx1].size(); idx2++) {
      size_t dummy = 0;
      m_process_S(idx1, idx2, dummy, true);
    }
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_DEC<T>::m_process_S(size_t idx1, size_t idx2, size_t& counter, bool read)
{
  auto& set = m_LIS[idx1][idx2];
  bool is_sig = true;

  if (read)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;  // Let's increment the counter first!
    m_code_S(idx1, idx2);
    set.set_length(0);  // this current set is gonna be discarded.
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_DEC<T>::m_process_P(size_t idx, size_t& counter, bool read)
{
  bool is_sig = true;

  if (read)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;  // Let's increment the counter first!
    m_sign_array.wbit(idx, m_bit_buffer.rbit());

    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

template <typename T>
inline void SZ3::SPERR::SPECK1D_INT_DEC<T>::m_code_S(size_t idx1, size_t idx2)
{
  auto subsets = m_partition_set(m_LIS[idx1][idx2]);
  auto sig_counter = size_t{0};
  auto read = bool{true};

  // Process the 1st subset
  const auto& set0 = subsets[0];
  assert(set0.get_length() != 0);
  if (set0.get_length() == 1) {
    m_LIP_mask.wtrue(set0.get_start());
    m_process_P(set0.get_start(), sig_counter, read);
  }
  else {
    const auto newidx1 = set0.get_level();
    m_LIS[newidx1].push_back(set0);
    m_process_S(newidx1, m_LIS[newidx1].size() - 1, sig_counter, read);
  }

  // Process the 2nd subset
  if (sig_counter == 0)
    read = false;
  const auto& set1 = subsets[1];
  assert(set1.get_length() != 0);
  if (set1.get_length() == 1) {
    m_LIP_mask.wtrue(set1.get_start());
    m_process_P(set1.get_start(), sig_counter, read);
  }
  else {
    const auto newidx1 = set1.get_level();
    m_LIS[newidx1].push_back(set1);
    m_process_S(newidx1, m_LIS[newidx1].size() - 1, sig_counter, read);
  }
}

// ---- END SPECK1D_INT_DEC.cpp ----

// ---- BEGIN Outlier_Coder.cpp ----
#include "Outlier_Coder.h"

#include <algorithm>
#include <cassert>
#include <cfenv>
#include <cfloat>  // FLT_ROUNDS
#include <cmath>

SZ3::SPERR::Outlier::Outlier(size_t p, double e) : pos(p), err(e) {}

inline void SZ3::SPERR::Outlier_Coder::add_outlier(Outlier out)
{
  m_LOS.push_back(out);
}

inline void SZ3::SPERR::Outlier_Coder::use_outlier_list(std::vector<Outlier> los)
{
  m_LOS = std::move(los);
}

inline void SZ3::SPERR::Outlier_Coder::set_length(size_t len)
{
  m_total_len = len;
}

inline void SZ3::SPERR::Outlier_Coder::set_tolerance(double tol)
{
  m_tol = tol;
}

inline auto SZ3::SPERR::Outlier_Coder::view_outlier_list() const -> const std::vector<Outlier>&
{
  return m_LOS;
}

inline void SZ3::SPERR::Outlier_Coder::append_encoded_bitstream(vec8_type& buf) const
{
  // Just append the bitstream produced by `m_encoder` is fine.
  std::visit([&buf](auto&& enc) { enc.append_encoded_bitstream(buf); }, m_encoder);
}

inline auto SZ3::SPERR::Outlier_Coder::get_stream_full_len(const void* p) const -> size_t
{
  return std::visit([p](auto&& dec) { return dec.get_stream_full_len(p); }, m_decoder);
}

inline auto SZ3::SPERR::Outlier_Coder::use_bitstream(const void* p, size_t len) -> RTNType
{
  // Decide on the integer length to use.
  const auto num_bitplanes = speck_int_get_num_bitplanes(p);
  if (num_bitplanes <= 8)
    m_instantiate_uvec_coders(UINTType::UINT8);
  else if (num_bitplanes <= 16)
    m_instantiate_uvec_coders(UINTType::UINT16);
  else if (num_bitplanes <= 32)
    m_instantiate_uvec_coders(UINTType::UINT32);
  else
    m_instantiate_uvec_coders(UINTType::UINT64);

  // Clean up data structures.
  m_sign_array.resize(0);
  m_LOS.clear();
  std::visit([](auto&& vec) { vec.clear(); }, m_vals_ui);

  // Ask the decoder to use the bitstream directly.
  std::visit([p, len](auto&& dec) { dec.use_bitstream(p, len); }, m_decoder);

  return RTNType::Good;
}

inline auto SZ3::SPERR::Outlier_Coder::encode() -> RTNType
{
  // Sanity check: whether we can proceed with encoding.
  if (m_total_len == 0 || m_tol <= 0.0 || m_LOS.empty())
    return RTNType::Error;
  if (std::any_of(m_LOS.cbegin(), m_LOS.cend(), [len = m_total_len, tol = m_tol](auto out) {
        return out.pos >= len || std::abs(out.err) <= tol;
      }))
    return RTNType::Error;

  // Step 1: find the biggest magnitude of outlier errors, and then instantiate data structures.
  auto maxerr = *std::max_element(m_LOS.cbegin(), m_LOS.cend(), [](auto v1, auto v2) {
    return std::abs(v1.err) < std::abs(v2.err);
  });
  std::fesetround(FE_TONEAREST);
  assert(FE_TONEAREST == std::fegetround());
  assert(FLT_ROUNDS == 1);
  std::feclearexcept(FE_INVALID);
  auto maxint = std::llrint(std::abs(maxerr.err));
  if (std::fetestexcept(FE_INVALID))
    return RTNType::FE_Invalid;

  if (maxint <= std::numeric_limits<uint8_t>::max())
    m_instantiate_uvec_coders(UINTType::UINT8);
  else if (maxint <= std::numeric_limits<uint16_t>::max())
    m_instantiate_uvec_coders(UINTType::UINT16);
  else if (maxint <= std::numeric_limits<uint32_t>::max())
    m_instantiate_uvec_coders(UINTType::UINT32);
  else
    m_instantiate_uvec_coders(UINTType::UINT64);

  // Step 2: quantize the outliers.
  m_quantize();

  // Step 3: integer SPECK encoding.
  std::visit([len = m_total_len](auto&& enc) { enc.set_dims({len, 1, 1}); }, m_encoder);
  auto rtn = RTNType::Good;
  switch (m_encoder.index()) {
    case 0:
      rtn = std::get<0>(m_encoder).use_coeffs(std::move(std::get<0>(m_vals_ui)),
                                              std::move(m_sign_array));
      break;
    case 1:
      rtn = std::get<1>(m_encoder).use_coeffs(std::move(std::get<1>(m_vals_ui)),
                                              std::move(m_sign_array));
      break;
    case 2:
      rtn = std::get<2>(m_encoder).use_coeffs(std::move(std::get<2>(m_vals_ui)),
                                              std::move(m_sign_array));
      break;
    default:
      rtn = std::get<3>(m_encoder).use_coeffs(std::move(std::get<3>(m_vals_ui)),
                                              std::move(m_sign_array));
  }
  if (rtn != RTNType::Good)
    return rtn;

  std::visit([](auto&& enc) { enc.encode(); }, m_encoder);

  return RTNType::Good;
}

inline auto SZ3::SPERR::Outlier_Coder::decode() -> RTNType
{
  // Sanity check
  if (m_total_len == 0 || m_tol <= 0.0)
    return RTNType::Error;

  // Step 1: Integer SPECK decode
  std::visit([len = m_total_len](auto&& dec) { dec.set_dims({len, 1, 1}); }, m_decoder);
  std::visit([](auto&& dec) { dec.decode(); }, m_decoder);
  std::visit([&vec = m_vals_ui](auto&& dec) { vec = dec.release_coeffs(); }, m_decoder);
  m_sign_array = std::visit([](auto&& dec) { return dec.release_signs(); }, m_decoder);

  // Step 2: Inverse quantization
  m_inverse_quantize();

  return RTNType::Good;
}

inline void SZ3::SPERR::Outlier_Coder::m_instantiate_uvec_coders(UINTType type)
{
  switch (type) {
    case UINTType::UINT8:
      if (m_vals_ui.index() != 0)
        m_vals_ui.emplace<0>();
      if (m_encoder.index() != 0)
        m_encoder.emplace<0>();
      if (m_decoder.index() != 0)
        m_decoder.emplace<0>();
      break;
    case UINTType::UINT16:
      if (m_vals_ui.index() != 1)
        m_vals_ui.emplace<1>();
      if (m_encoder.index() != 1)
        m_encoder.emplace<1>();
      if (m_decoder.index() != 1)
        m_decoder.emplace<1>();
      break;
    case UINTType::UINT32:
      if (m_vals_ui.index() != 2)
        m_vals_ui.emplace<2>();
      if (m_encoder.index() != 2)
        m_encoder.emplace<2>();
      if (m_decoder.index() != 2)
        m_decoder.emplace<2>();
      break;
    default:
      if (m_vals_ui.index() != 3)
        m_vals_ui.emplace<3>();
      if (m_encoder.index() != 3)
        m_encoder.emplace<3>();
      if (m_decoder.index() != 3)
        m_decoder.emplace<3>();
  }
}

inline void SZ3::SPERR::Outlier_Coder::m_quantize()
{
  std::visit([len = m_total_len](auto&& vec) { vec.assign(len, 0); }, m_vals_ui);
  m_sign_array.resize(m_total_len);
  m_sign_array.reset_true();

  std::visit(
      [&los = m_LOS, &signs = m_sign_array, tol = m_tol](auto&& vec) {
        auto inv = 1.0 / tol;
        for (auto out : los) {
          auto ll = std::llrint(out.err * inv);
          signs.wbit(out.pos, ll >= 0);
          vec[out.pos] = std::abs(ll);
        }
      },
      m_vals_ui);
}

inline void SZ3::SPERR::Outlier_Coder::m_inverse_quantize()
{
  m_LOS.clear();

  // First, bring all non-zero integer correctors to `m_LOS`.
  std::visit(
      [&los = m_LOS](auto&& vec) {
        for (size_t i = 0; i < vec.size(); i++)
          switch (vec[i]) {
            case 0:
              break;
            case 1:
              los.emplace_back(i, 1.1);
              break;
            default:
              los.emplace_back(i, static_cast<double>(vec[i]) - 0.25);
          }
      },
      m_vals_ui);

  // Second, restore the floating-point correctors.
  const auto tmp = std::array<double, 2>{-1.0, 1.0};
  std::transform(m_LOS.cbegin(), m_LOS.cend(), m_LOS.begin(),
                 [q = m_tol, &signs = m_sign_array, tmp](auto los) {
                   auto b = signs.rbit(los.pos);
                   los.err *= (q * tmp[b]);
                   return los;
                 });
}
// ---- END Outlier_Coder.cpp ----

#endif
