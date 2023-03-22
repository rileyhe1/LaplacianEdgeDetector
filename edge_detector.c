#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <string.h>
//change the number of threads as you run your concurrency experiment
#define LAPLACIAN_THREADS 500

/* Laplacian filter is 3 by 3 */
#define FILTER_WIDTH 3
#define FILTER_HEIGHT 3

#define RGB_COMPONENT_COLOR 255
#define MAX_RGB_COLOR 255
#define MIN_RGB_COLOR 0
#define SIZE_OF_MAX_COLOR_VALUE_STR 4

#define ONE_ARG 2
#define BUFFER_SIZE 100
#define SIZE_OF_IMG_FORMAT_STR 3
#define MAX_FILENAME_LENGTH 100
#define MICRO_TO_MILLI 0.001
#define MILLI_TO_SECONDS 0.001
#define ONE_WRITE 1
#define ONE_READ 1
#define SIZE_OF_ONE_CHAR 1
#define ERR_NUM 1
#define SUCCESS 0

typedef struct
{
      unsigned char r, g, b;
} PPMPixel;

struct parameter
{
    PPMPixel *image;         //original image pixel data
    PPMPixel *result;        //filtered image pixel data
    unsigned long int w;     //width of image
    unsigned long int h;     //height of image
    unsigned long int start; //starting point of work
    unsigned long int size;  //equal share of work (almost equal if odd)
};


struct file_name_args
{
    char *input_file_name;      //e.g., file1.ppm
    char output_file_name[20];  //will take the form laplaciani.ppm, e.g., laplacian1.ppm
};

double total_elapsed_time = 0;


/*  This is the thread function. It will compute the new values for the region of image
    specified in params (start to start+size) using convolution. For each pixel in the
    input image, the filter is conceptually placed on top of the image with its origin
    lying on that pixel. The  values  of  each  input  image  pixel  under  the  mask
    are  multiplied  by the corresponding filter values. Truncate values smaller than
    zero to zero and larger than 255 to 255. The results are summed together to yield
    a single output value that is placed in the output image at the location of the
    pixel being processed on the input.

    First compute_laplacian_threadfn casts params back to a parameter pointer, and then
    iterates over each pixel in its territory and does the laplacian math, and returns the result.
 */
void *compute_laplacian_threadfn(void *params)
{
    int laplacian[FILTER_WIDTH][FILTER_HEIGHT] =
    {
        {-1, -1, -1},
        {-1,  8, -1},
        {-1, -1, -1}
    };

    int red, green, blue;
    struct parameter* p = (struct parameter*)params;
    for(int iteratorImageHeight = p->start; iteratorImageHeight < p->size + p->start; iteratorImageHeight++)
    {
      for(int iteratorImageWidth = 0; iteratorImageWidth < p->w; iteratorImageWidth++)
      {
        int x_coordinate = 0;
        int y_coordinate = 0;
        red = 0;
        blue = 0;
        green = 0;
        for(int iteratorFilterHeight = 0; iteratorFilterHeight < FILTER_HEIGHT; iteratorFilterHeight++)
        {
          for(int iteratorFilterWidth = 0; iteratorFilterWidth < FILTER_WIDTH; iteratorFilterWidth++)
          {
            x_coordinate = (iteratorImageWidth - FILTER_WIDTH / 2 + iteratorFilterWidth + p->w) % p->w;
            y_coordinate = (iteratorImageHeight - FILTER_HEIGHT / 2 + iteratorFilterHeight + p->h) % p->h;
            red += p->image[y_coordinate * p->w + x_coordinate].r * laplacian[iteratorFilterHeight][iteratorFilterWidth];
            green+= p->image[y_coordinate * p->w + x_coordinate].g * laplacian[iteratorFilterHeight][iteratorFilterWidth];
            blue+= p->image[y_coordinate * p->w + x_coordinate].b * laplacian[iteratorFilterHeight][iteratorFilterWidth];
          }
        }
        if(red > MAX_RGB_COLOR)
        {
          red = MAX_RGB_COLOR;
        }
        else if (red < MIN_RGB_COLOR)
        {
          red = MIN_RGB_COLOR;
        }
        if(green > MAX_RGB_COLOR)
        {
          green = MAX_RGB_COLOR;
        }
        else if (green < MIN_RGB_COLOR)
        {
          green = MIN_RGB_COLOR;
        }
        if(blue > MAX_RGB_COLOR)
        {
          blue = MAX_RGB_COLOR;
        }
        else if (blue < MIN_RGB_COLOR)
        {
          blue = MIN_RGB_COLOR;
        }
        p->result[iteratorImageHeight * p->w + iteratorImageWidth].r = red;
        p->result[iteratorImageHeight * p->w + iteratorImageWidth].g = green;
        p->result[iteratorImageHeight * p->w + iteratorImageWidth].b = blue;
      }
    }
    return NULL;
}

/* Apply the Laplacian filter to an image using threads.
 Each thread shall do an equal share of the work, i.e. work=height/number of threads.
 If the size is not even, the last thread shall take the rest of the work.
 Compute the elapsed time and store it in *elapsedTime (Read about gettimeofday).
 Return: result (filtered image)

    apply_filters first allocates space for the result image, and then for each thread's
    parameter. Next it calculates how much work each thread should be doing, with any leftover
    work done by the final thread. It sets the values for each of the parameter's members, and
    finally creates each thread, joins them and frees the resources it allocated at the end.
 */
PPMPixel *apply_filters(PPMPixel *image, unsigned long w, unsigned long h, double *elapsedTime)
{

    PPMPixel *result = (PPMPixel*)malloc(w * h * sizeof(PPMPixel));
    struct parameter* p = (struct parameter*)malloc(LAPLACIAN_THREADS * sizeof(struct parameter));
    int work = h / LAPLACIAN_THREADS;
    int remainder = h % LAPLACIAN_THREADS;

    pthread_t threads[LAPLACIAN_THREADS];
    for(int i = 0; i < LAPLACIAN_THREADS; i++)
    {
      p->image = image;
      p->result = result;
      p->w = w;
      p->h = h;
      p->start = (unsigned long int) 0 + (work * i);
      if(i == LAPLACIAN_THREADS -1)
      {
        p->size = (unsigned long int)work + remainder;
      }
      else
      {
        p->size = (unsigned long int) work;
      }
      pthread_create(&threads[i], NULL, compute_laplacian_threadfn, (void*)p);
      p++;
    }
    for(int i = 0; i < LAPLACIAN_THREADS; i++)
    {
      pthread_join(threads[i], NULL);
      p--;
    }
    free(p);
    return result;
}

/*Create a new P6 file to save the filtered image in. Write the header block
 e.g. P6
      Width Height
      Max color value
 then write the image data.
 The name of the new file shall be "filename" (the second argument).

    write_image tries to open the output file and prints an error if it fails, then writes the image
    header information followed by the pixel information.
 */
void write_image(PPMPixel *image, char *filename, unsigned long int width, unsigned long int height)
{
  FILE* ofp  = fopen(filename, "w+");
  if(ofp == NULL)
  {
    perror("error opening output file...");
    return;
  }

  fwrite("P6\n", SIZE_OF_IMG_FORMAT_STR, ONE_WRITE, ofp);
  char intStrW [BUFFER_SIZE];
  char intStrH [BUFFER_SIZE];
  sprintf(intStrW, "%ld", width);
  sprintf(intStrH, "%ld", height);
  fwrite(intStrW, strlen(intStrW), ONE_WRITE, ofp);
  fwrite(" ", SIZE_OF_ONE_CHAR, ONE_WRITE, ofp);
  fwrite(intStrH, strlen(intStrH), ONE_WRITE, ofp);
  fwrite("\n", SIZE_OF_ONE_CHAR, ONE_WRITE, ofp);
  fwrite("255\n", SIZE_OF_MAX_COLOR_VALUE_STR, ONE_WRITE, ofp);

  int numPixels = width*height;
  for(int i = 0; i < numPixels; i++)
  {
    fwrite(image, sizeof(PPMPixel), ONE_WRITE, ofp);
    image++;
  }
  fclose(ofp);
}



/* Open the filename image for reading, and parse it.
    Example of a ppm header:    //http://netpbm.sourceforge.net/doc/ppm.html
    P6                  -- image format
    # comment           -- comment lines begin with
    ## another comment  -- any number of comment lines
    200 300             -- image width & height
    255                 -- max color value

 Check if the image format is P6. If not, print invalid format error message.
 If there are comments in the file, skip them. You may assume that comments exist
 only in the header block. Read the image size information and store them in width
 and height. Check the rgb component, if not 255, display error message.
 Return: pointer to PPMPixel that has the pixel data of the input image (filename).
 The pixel data is stored in scanline order from left to right (up to bottom) in 3-byte
 chunks (r g b values for each pixel) encoded as binary numbers.
 */
 /*
 This function does a series of error checking to ensure the file is a ppm file, in P6 format,
 has a maximum color value of 255, and checks the image file's width and height and stores those values
 in width and height respecively. If the file does not satisfy those requirements, an error message is printed
 and NULL is returned. Next the function reads the pixel data in its entirety and stores it in the img pointer,
 after mallocing appropriate space for it based on the width and height of the image.
 */
PPMPixel *read_image(const char *filename, unsigned long int *width, unsigned long int *height)
{
    char filetype [BUFFER_SIZE];
    int locationOfDot = strcspn(filename, ".");
    int fileTypeIdx = 0;

    for(int i = locationOfDot+1; i < strlen(filename)+1; i++)
    {
      filetype[fileTypeIdx] = filename[i];
      fileTypeIdx++;
    }
    if(strcmp("ppm", filetype) != 0)
    {
      perror("Error: input must be ppm file\n");
      return NULL;
    }

    FILE* fp = fopen(filename, "rb");
    if(fp == NULL)
    {
      perror("Error opening input file...\n");
      return NULL;
    }

    char imgFormat [BUFFER_SIZE];
    if(fgets(imgFormat, BUFFER_SIZE, fp) == NULL)
    {
      perror("error reading input file\n");
      return NULL;
    }
    int locationOfNewline = strcspn(imgFormat, "\n");
    imgFormat[locationOfNewline] = '\0';
    if(strcmp(imgFormat, "P6") != 0)
    {
      perror("Image Format must be P6...\n");
      return NULL;
    }


    char imgDimensions[BUFFER_SIZE];
    fgets(imgDimensions, BUFFER_SIZE, fp);
    while(imgDimensions[0] == '#')
    {
      fgets(imgDimensions, BUFFER_SIZE, fp);
    }
    locationOfNewline = strcspn(imgDimensions, "\n");
    imgDimensions[locationOfNewline] = '\0';
    int locationOfSpace = strcspn(imgDimensions, " ");
    char widthStr[BUFFER_SIZE];
    char heightStr[BUFFER_SIZE];
    int widthIdx = 0;
    int heightIdx = 0;
    for(int i = 0; i < strlen(imgDimensions)+1; i++)
    {
      if(i < locationOfSpace)
      {
        widthStr[widthIdx] = imgDimensions[i];
        widthIdx++;
      }
      else if(i > locationOfSpace)
      {
        heightStr[heightIdx] = imgDimensions[i];
        heightIdx++;
      }
    }
    widthStr[widthIdx] = '\0';
    *width = (unsigned long int) atoi(widthStr);
    *height = (unsigned long int) atoi(heightStr);

    char maxColorValue[BUFFER_SIZE];
    fgets(maxColorValue, BUFFER_SIZE, fp);
    locationOfNewline = strcspn(maxColorValue, "\n");
    maxColorValue[locationOfNewline] = '\0';
    if(atoi(maxColorValue) != RGB_COMPONENT_COLOR)
    {
      perror("Invalid maximum color value for input file...\n");
      return NULL;
    }

    PPMPixel *img = (PPMPixel*)malloc((*width)*(*height) * sizeof(PPMPixel));
    int numPixels = (*width) * (*height);

    for(int i = 0; i < numPixels; i++)
    {
      fread(img, sizeof(PPMPixel), ONE_READ, fp);
      img++;
    }

    for(int i = 0; i < numPixels; i++)
    {
      img--;
    }

    fclose(fp);
    return img;
}

/* The thread function that manages an image file.
 Read an image file that is passed as an argument at runtime.
 Apply the Laplacian filter.
 Save the result image in a file called laplaciani.ppm, where i is the image file order
in the passed arguments. Example: the result image of the file passed third during the
input shall be called "laplacian3.ppm".

    manage_image_file first casts args back into a file_name_args pointer, and mallocs a pointer
    for both width and height of the input image, and calls read, apply_filters, and write_image,
    and lastly frees the resources it allocated at the start.
    If input image is NULL, this means something went wrong in the read function, so we want the thread
    to exit.
*/
void *manage_image_file(void *args)
{
  struct file_name_args* arguments = (struct file_name_args*)args;
  unsigned long int* width = (unsigned long int*) malloc(sizeof(unsigned long int*));
  unsigned long int* height = (unsigned long int*) malloc(sizeof(unsigned long int*));
  PPMPixel* inputImg;
  PPMPixel* outputImg;

  // read
  inputImg = read_image(arguments->input_file_name, width, height);

  if(inputImg == NULL)
  {
    exit(ERR_NUM);
  }
  // apply filters
  outputImg = apply_filters(inputImg, *width, *height, &total_elapsed_time);
  // write
  write_image(outputImg, arguments->output_file_name, *width, *height);
  free(inputImg);
  free(outputImg);
  free(width);
  free(height);

  return NULL;
}
/*The driver of the program. Check for the correct number of arguments. If wrong print
 the message: "Usage ./a.out filename[s]" It shall accept n filenames as arguments,
 separated by whitespace, e.g., ./a.out file1.ppm file2.ppm    file3.ppm
  It will create a thread for each input file to manage.
  It will print the total elapsed time in .4 precision seconds(e.g., 0.1234 s).
  The total elapsed time is the total time taken by all threads to compute the edge detection
  of all input images .
 */

 /*
    main initializes an array of threads and mallocs the file_name_args pointers for each respective
    thread. it then creates a string to hold i as a string which it converts using sprintf. After
    each thread is created it uses pthread_join to wait for all the threads to finish executing.
    Lastly it frees the allocated file_name_args pointers and prints the total elapsed time in seconds.
 */
  int main(int argc, char *argv[])
  {
      if(argc < ONE_ARG)
      {
        printf("Usage ./a.out filename\n");
        return -1;
      }

      struct timeval timer1, timer2;
      gettimeofday(&timer1, NULL);

      pthread_t threads[argc -1];
      struct file_name_args* args = (struct file_name_args*)malloc((argc-1) * sizeof(struct file_name_args));

      for(int i = 0; i < argc-1; i++)
      {
        char filename[MAX_FILENAME_LENGTH];
        char iAsStr[MAX_FILENAME_LENGTH];
        strcpy(filename, "laplacian");
        sprintf(iAsStr, "%d", i+1);
        strcat(filename, iAsStr);
        strcat(filename, ".ppm");
        args->input_file_name = argv[i+1];
        strcpy(args->output_file_name, filename);
        pthread_create(&threads[i], NULL, manage_image_file, (void*)args);
        args++;
      }
      for(int i = 0; i < argc-1; i++)
      {
        pthread_join(threads[i], NULL);
        args--;
      }
      free(args);
      gettimeofday(&timer2, NULL);
      total_elapsed_time += (timer2.tv_sec - timer1.tv_sec);
      total_elapsed_time += ((timer2.tv_usec - timer1.tv_usec)* MICRO_TO_MILLI) * MILLI_TO_SECONDS;
      printf("Total elapsed time: %.3f\n", total_elapsed_time);
      return SUCCESS;
  }
