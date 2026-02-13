#!/bin/bash

doxygen docs/Doxyfile
npx -y pagefind --site docs/html
