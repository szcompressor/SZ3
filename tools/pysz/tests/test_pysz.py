"""
Test script for pysz - Python interface for SZ3 compression

Usage:
    pytest tests/                    # Run with pytest
    python tests/test_pysz.py       # Run standalone
"""

import sys
import numpy as np

from pysz import sz, pyConfig


def test_compression():
    """Test basic compression and decompression"""
    
    print("=" * 70)
    print("pysz Basic Test")
    print("=" * 70)
    
    # Create test data
    print("\n[1] Creating test data...")
    data = np.random.randn(100, 100).astype(np.float32)
    print(f"    ✓ Shape: {data.shape}, dtype: {data.dtype}, size: {data.nbytes} bytes")
    
    # Create config
    print("[2] Creating config...")
    config = pyConfig(100, 100)
    config.errorBoundMode = pyConfig.EB.ABS
    config.absErrorBound = 0.01
    print(f"    ✓ Mode: {config.errorBoundMode}, Bound: {config.absErrorBound}")
    
    # Compress
    print("[3] Compressing...")
    compressed, ratio = sz.compress(data, config)
    print(f"    ✓ Ratio: {ratio:.2f}x ({data.nbytes} → {len(compressed)} bytes)")
    
    # Decompress
    print("[4] Decompressing...")
    decompressed = sz.decompress(compressed, data.dtype, data.shape, config)
    print(f"    ✓ Shape: {decompressed.shape}")
    
    # Verify
    print("[5] Verifying...")
    max_error, psnr, nrmse = sz.verify(data, decompressed)
    print(f"    ✓ Max error: {max_error:.2e}, PSNR: {psnr:.2f} dB, NRMSE: {nrmse:.2e}")
    
    assert max_error <= 0.01, f"Error {max_error} exceeds bound 0.01"
    
    # Test double precision
    print("[6] Testing double precision...")
    data_double = np.random.randn(50, 50).astype(np.float64)
    config_double = pyConfig(50, 50)
    config_double.errorBoundMode = pyConfig.EB.ABS
    config_double.absErrorBound = 1e-6
    compressed_d, ratio_d = sz.compress(data_double, config_double)
    decompressed_d = sz.decompress(compressed_d, data_double.dtype, data_double.shape, config_double)
    max_error_d, _, _ = sz.verify(data_double, decompressed_d)
    print(f"    ✓ Ratio: {ratio_d:.2f}x, Max error: {max_error_d:.2e}")
    assert max_error_d <= 1e-6, f"Double error {max_error_d} exceeds bound 1e-6"
    
    # Test 3D data
    print("[7] Testing 3D data...")
    data_3d = np.random.randn(20, 30, 40).astype(np.float32)
    config_3d = pyConfig(20, 30, 40)
    config_3d.errorBoundMode = pyConfig.EB.REL
    config_3d.relErrorBound = 0.001
    compressed_3d, ratio_3d = sz.compress(data_3d, config_3d)
    decompressed_3d = sz.decompress(compressed_3d, data_3d.dtype, data_3d.shape, config_3d)
    max_error_3d, _, _ = sz.verify(data_3d, decompressed_3d)
    print(f"    ✓ Ratio: {ratio_3d:.2f}x, Max error: {max_error_3d:.2e}")
    
    print("\n" + "=" * 70)
    print("All tests passed!")
    print("=" * 70)


if __name__ == '__main__':
    try:
        test_compression()
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
