#include <nemo.h>
#include <gnuastro/fits.h>


string defv[] = {
    "in=???\n           Input fits file",
    "hdu=0\n            The HDU number to pick (0 is the first)",
    "VERSION=0.2\n      25-oct-2021 PJT",
    NULL,
};

string usage = "read fits file using gnuastro (demo1)";


// taken from https://www.gnu.org/software/gnuastro/manual/html_node/Library-demo-_002d-reading-a-image.html
void nemo_main(void)
{
  size_t i;
  float *farray;
  double sum=0.0f;
  gal_data_t *image;
  char *filename = getparam("in");
  char *hdu = getparam("hdu");          // counter intuitive, why is not this an int?

  /* Read `img.fits' (HDU: 1) as a float32 array. */
  image=gal_fits_img_read_to_type(filename, hdu, GAL_TYPE_FLOAT32, -1, 1);

  /* Use the allocated space as a single precision floating
   * point array (recall that `image->array' has `void *'
   * type, so it is not directly usable. */
  farray=image->array;


  /* Calculate the sum of all the values. */
  for(i=0; i<image->size; ++i)
    sum += farray[i];

  /* Report the sum. */
  printf("Sum of values in %s (hdu %s) is: %f\n", filename, hdu, sum);

  /* Clean up and return. */
  gal_data_free(image);
  // return EXIT_SUCCESS;
}
