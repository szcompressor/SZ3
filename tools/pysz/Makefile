.PHONY: install build clean test help

help:
	@echo "pysz - Python bindings for SZ3"
	@echo ""
	@echo "Available targets:"
	@echo "  install    - Install package in editable mode"
	@echo "  build      - Build the package"
	@echo "  clean      - Remove all build artifacts"
	@echo "  test       - Run tests"
	@echo "  help       - Show this help message"

install:
	pip install -e .

build:
	python -m build

clean:
	@echo "Cleaning build artifacts..."
	rm -rf build/ dist/ .eggs/ *.egg-info
	rm -rf src/*.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*.so" -delete
	find . -type f -name "*.dylib" -delete
	find . -type f -name "*.cpp" -path "*/pysz/*" -delete
	find . -type f -name "*.c" -path "*/pysz/*" -delete
	find . -type f -name "*.html" -delete
	@echo "✓ Cleaned!"

test:
	python tests/test_basic.py
