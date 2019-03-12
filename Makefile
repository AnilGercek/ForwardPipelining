rasterizer_cpp:
	g++ -o rasterizer hw2_math_ops.cpp hw2_file_ops.cpp rasterizer.cpp -g

clean:
	rm -rf *.ppm *.png ./rasterizer
