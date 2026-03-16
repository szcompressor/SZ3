with import <nixpkgs> { };
stdenv.mkDerivation {
  pname = "sz3";
  version = "3.3.x";

  src = fetchFromGitHub {
    owner = "szcompressor";
    repo = "SZ3";
    rev = "34b36f2af0d43aa4f5dea0e39500e91f6d35c9cf";
    sha256 = "0fckf5939m9d2w3095isnq3qfrz8fxpppn7y3c2hwzdh72vdzic7";
  };

  nativeBuildInputs = [
    cmake
    pkg-config
  ];

  buildInputs = [
    zstd
    gsl
  ]
  ++ lib.optionals stdenv.cc.isClang [
    llvmPackages.openmp
  ];

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
  ];

  meta = with lib; {
    description = "A modular error-bounded lossy compression framework for scientific datasets";
    homepage = "https://github.com/szcompressor/SZ3";
    maintainer = [ "szcomprssor" ];
    license = licenses.bsd3;
    platforms = platforms.unix;
  };
}
