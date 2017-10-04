# Kernel-Based-Image-Processing

The program implements three kernel based image processing filters, blurr, sobel and sharp.


To compile, simply type "make" <br>
To run, ./filter "input bmp fileName" "filter parameters" "output bmp fileName" <br>

Example, for blur ./filter tommy.bmp blur 5 2.0 tommy_blur.bmp
<br> for sobel
./filter bw_tile_wikimedia.bmp sobel bw_tile_sobel.bmp
<br> for sharp
./filter usc_ucla_wikimedia.bmp unsharp 5 2.0 0.7 usc_ucla_unsharp.bmp
