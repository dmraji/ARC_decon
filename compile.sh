echo "Removing old build"
rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt decon_main Makefile

echo "Building..."
cmake ~/4232/radiation_processing

echo "Done! Now run $ make to compile."
