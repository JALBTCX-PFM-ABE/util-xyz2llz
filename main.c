
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "nvutility.h"

#include "llz.h"

#include "version.h"


/*****************************************************************************

    Program:    xyz2llz

    Purpose:    Converts PFM bin or index data to llz format for input to the
                chrtr gridding program.

    Programmer: Jan C. Depner

    Date:       11/08/99

*****************************************************************************/


void usage ()
{
  fprintf (stderr, "\nUsage: xyz2llz [-i] [-t] INPUT_FILENAME OUTPUT_FILENAME\n");
  fprintf (stderr, "\nWhere:\n\n");
  fprintf (stderr, "\tINPUT_FILENAME must end in .xyz or .yxz depending on file type.\n");
  fprintf (stderr, "\t-i  =  invert sign of input Z value.\n");
  fprintf (stderr, "\t-t  =  input file has time as the 4th field (after Z).\n\n");
  fprintf (stderr, "\n\tInput data for .yxz must be lat,lon,z or lat lon z.\n");
  fprintf (stderr, "\n\tFor .xyz it must be lon,lat,z or lon lat z (all positions in degrees).\n");
  fprintf (stderr, "\n\tIMPORTANT NOTE: Time must be in decimal seconds from 01/01/1970.\n\n");
  fflush (stderr);
}



int32_t main (int32_t argc, char **argv)
{
  FILE                *fp = NULL;
  int32_t             percent = 0, old_percent = -1, out_count = 0, llz_hnd = -1, size = 0;
  float               z;
  char                out_file[256], in_file[256], string[512];
  uint8_t             yxz = NVFalse, flip = NVFalse, tflag = NVFalse;
  char                c;
  double              lat, lon, t;
  LLZ_REC             llz_rec;
  LLZ_HEADER          llz_header;
  extern char         *optarg;
  extern int          optind;



  printf ("\n\n %s \n\n\n", VERSION);

  yxz = NVFalse;
  flip = NVFalse;
  tflag = NVFalse;

  while ((c = getopt (argc, argv, "it")) != EOF)
    {
      switch (c)
        {
        case 'i':
          flip = NVTrue;
          break;

        case 't':
          tflag = NVTrue;
          break;

        default:
          usage ();
          exit (-1);
          break;
        }
    }


  /* Make sure we got the mandatory file name arguments.  */

  if (optind >= argc || (!strstr (argv[optind], ".xyz") && !strstr (argv[optind], ".yxz")))
    {
      usage ();
      exit (-1);
    }


  if (strstr (argv[optind], ".yxz")) yxz = NVTrue;


  strcpy (in_file, argv[optind]);


  if ((fp = fopen (in_file, "r")) == NULL)
    {
      perror (in_file);
      exit (-1);
    }


  fprintf (stderr, "File : %s\n\n", in_file);
  fflush (stderr);


  /*  Figure out the file size.  */

  fseek (fp, 0, SEEK_END);
  size = ftell (fp);
  fseek (fp, 0, SEEK_SET);


  strcpy (out_file, argv[optind + 1]);

  if (strcmp (&out_file[strlen (out_file) - 4], ".llz")) strcat (out_file, ".llz");


  /*  Boilerplate LLZ header.  */

  sprintf (llz_header.comments, "Created from %s using %s", gen_basename (in_file), VERSION);
  llz_header.time_flag = tflag;
  llz_header.depth_units = 0;

  if ((llz_hnd = create_llz (out_file, llz_header)) < 0)
    {
      perror (out_file);
      exit (-1);
    }


  out_count = 0;
  while (ngets (string, 512, fp) != NULL)
    {
      /*  Skip comments.  */

      if (string[0] != '#')
        {
          llz_rec.tv_sec = 0;
          llz_rec.tv_nsec = 0;


          if (yxz)
            {
              if (strchr (string, ','))
                {
                  if (tflag)
                    {
                      sscanf (string, "%lf,%lf,%f,%lf", &lat, &lon, &z, &t);

                      llz_rec.tv_sec = (time_t) t;
                      llz_rec.tv_nsec = NINT (fmod (t, 1.0L) * 1000000000.0);
                    }
                  else
                    {
                      sscanf (string, "%lf,%lf,%f", &lat, &lon, &z);
                    }
                }
              else
                {
                  if (tflag)
                    {
                      sscanf (string, "%lf %lf %f %lf", &lat, &lon, &z, &t);

                      llz_rec.tv_sec = (time_t) t;
                      llz_rec.tv_nsec = NINT (fmod (t, 1.0L) * 1000000000.0);
                    }
                  else
                    {
                      sscanf (string, "%lf %lf %f", &lat, &lon, &z);
                    }
                }
            }
          else
            {
              if (strchr (string, ','))
                {
                  if (tflag)
                    {
                      sscanf (string, "%lf,%lf,%f,%lf", &lon, &lat, &z, &t);

                      llz_rec.tv_sec = (time_t) t;
                      llz_rec.tv_nsec = NINT (fmod (t, 1.0L) * 1000000000.0);
                    }
                  else
                    {
                      sscanf (string, "%lf,%lf,%f", &lon, &lat, &z);
                    }
                }
              else
                {
                  if (tflag)
                    {
                      sscanf (string, "%lf %lf %f %lf", &lon, &lat, &z, &t);

                      llz_rec.tv_sec = (time_t) t;
                      llz_rec.tv_nsec = NINT (fmod (t, 1.0L) * 1000000000.0);
                    }
                  else
                    {
                      sscanf (string, "%lf %lf %f", &lon, &lat, &z);
                    }
                }
            }

          if (flip) z = -z;

          llz_rec.xy.lat = lat;
          llz_rec.xy.lon = lon;
          llz_rec.depth = z;
          llz_rec.status = 0;

          append_llz (llz_hnd, llz_rec);

          out_count++;
        }


      percent = ((float) ftell (fp) / (float) size) * 100.0;
      if (percent != old_percent)
        {
          old_percent = percent;
          fprintf (stderr, "%03d%% processed            \r", percent);
          fflush (stderr);
        }
    }


  close_llz (llz_hnd);


  fprintf (stderr, "%d LLZ records output        \n\n\n", out_count);
  fflush (stderr);


  return (0);
}
