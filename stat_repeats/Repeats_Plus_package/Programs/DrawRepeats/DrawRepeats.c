/******************************************************************************
**
** Source File Name = draw_repeats.c
**
** Program call: draw_repeats input_dataset output_dataset left_interval right_interval [[g]el <sequence_length>] [s <sequence_string>]
**
**
** Program create LaTex dataset for drawing repeat sequences connection on x-axis
** Parameters: 
** left_interval right_interval - Compulsory. Represent "window" form left_interval and _right interval for drawing. Only 
**                                objects where both parts (left and right) lies inside this interval will be drawn
** [g]el <sequence_length>      - Optional. If exists, only sequence with length great or equal (if gel) or equal (el) to paramtere will be drawn
** s <sequence_string>          - Oprional. If exists, only sequences with left part equal to <sequence_string> will be drawn
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


FILE *input,*iztex;

int main(int argc, char *argv[])

{
 int  left_interval, right_interval=0;
 int  sequence_length =0;
 int  i; 
 char sequence[251];
 short gel=0;
 short el=0;
 short niska=0;
 char gelc[101];
 char lens[101];
 char ulaz[44];
 char izlaz[44];
 char line[801];
 char protein[51];
 int  levo1,desno1,levo2,desno2,duzina;
 float levo_skalirano,desni_skalirano, skalirana_vrednost1, skalirana_vrednost2;
 char leva_niska[251],desna_niska[251];
 int brojacboje=0;
 char* readLine;
char boja[21][24]={" ",
"black","blue","brown","cyan","darkgray",
"-cyan","green","-orange","lime","magenta",
"olive","orange","pink","purple","red",
"teal","violet","yellow","red!75!green!50!blue",
"red!75!green"
} ;
char aar[21][2]={" ","W","F","Y","I","M","L","V","N","C","T","A","G","R","D","H","Q","K","S","E","P"};

 if ((argc <5) || (argc >9) || (argc==6) || (argc==8))
  {
   printf("\n Wrong number of arguments. Number of arguments are in interval [4,8]. \n Program is invoked with command:\n draw_repeats input_dataset output_dataset left_interval right_interval [[g]el <sequence_length>] [s <sequence_string>]");
   printf("\n");
   return 1;
  }

 strcpy(ulaz,argv[1]);
 strcpy(izlaz,argv[2]);


 if ((input = fopen(ulaz,"r")) == NULL)
    {
     printf("\nDataset %s could not be open or note exists",ulaz);
     return 4;
    };    

 if ((iztex = fopen(izlaz,"w")) == NULL)
    {
     printf("\nDataset %s could not be open or note exists",izlaz);
     return 5;
    };    
 left_interval=atoi(argv[3]);
 right_interval=atoi(argv[4]);
 strcpy(sequence,"");
  
 if (argc>5)
    {
     if (argc==7)
        {
         strcpy(gelc,argv[5]);
         if (strcmp(gelc,"gel")==0)
            {
             gel=1;
             sequence_length=atoi(argv[6]);
            }
            else if (strcmp(gelc,"el")==0)
                    {
                     el=1;
                     sequence_length=atoi(argv[6]);
                    }
                    else if (strcmp(gelc,"s")==0)
                            {
                             niska=1;
                             strcpy(sequence,argv[6]);
                            }
                           else {
                                  printf("\nInvalid parameter %s\n",argv[5]);
                                  return 99;
                                }
          
        }
        else if (argc==9)
                {
                 strcpy(gelc,argv[5]);

                 if (strcmp(gelc,"gel")==0)
                    {
                     gel=1;
                     sequence_length=atoi(argv[6]);
                     strcpy(gelc,argv[7]);
                     if (strcmp(gelc,"s")==0)
                        {
                         niska=1;
                         strcpy(sequence,argv[8]);
                        }
                        else {
                              printf("\nInvalid parameter %s\n",argv[7]);
                              return 99;
                             } 
                    }
                    else if (strcmp(gelc,"el")==0)
                            {
                             el=1;
                             sequence_length=atoi(argv[6]);
                             strcpy(gelc,argv[7]);			     
                             if (strcmp(gelc,"s")==0)
                                {
                                 niska=1;
                                 strcpy(sequence,argv[8]);
                                }
                                else {
                                      printf("\nInvalid parameter %s\n",argv[7]);
                                      return 99;
                                     } 
                            }
                            else if (strcmp(gelc,"s")==0)
                                    {
                                     niska=1;
                                     strcpy(sequence,argv[6]);
                                     strcpy(gelc,argv[7]);
                                     if (strcmp(gelc,"gel")==0)
                                        {
                                         gel=1;
                                         sequence_length=atoi(argv[8]);
                                        }
                                        else if (strcmp(gelc,"el")==0)
                                                {
                                                 el=1;
                                                 sequence_length=atoi(argv[8]);
                                                }
                                                else {
                                                      printf("\nInvalid parameter %s\n",argv[7]);
                                                      return 99;
                                                     }                                                 
                                    
                                    }                                                                      
                                    else {
                                          printf("\nInvalid parameter %s\n",argv[7]);
                                          return 99;
                                         } 
                }
                else {
                      printf("\nInvalid number of parameters\n");
                      return 98;
                     }




    }
 fprintf(iztex,"\\documentclass[a4paper]{article}                                   ");
 fprintf(iztex,"\n\\usepackage{pgfplots}                                            ");
 fprintf(iztex,"\n\\usepackage{pdflscape}                                           ");
 fprintf(iztex,"\n\\hoffset -2.5cm                                                  ");
 fprintf(iztex,"\n\\voffset 14cm                                                    ");
 fprintf(iztex,"\n\\textwidth 260mm                                                 ");
 fprintf(iztex,"\n\\textheight 90mm                                                 ");
 fprintf(iztex,"\n\\pagestyle{empty}                                                ");
 fprintf(iztex,"\n\\pgfplotsset{width=26cm,compat=1.9}                              ");
 fprintf(iztex,"\n\\pagestyle{empty}                                                ");
 fprintf(iztex,"\n\\begin{document}                                                 ");
 fprintf(iztex,"\n\\thispagestyle{empty}                                            ");
 fprintf(iztex,"\n\\begin{landscape}                                                ");
 fprintf(iztex,"\n\\begin{tikzpicture}                                              "); 
 fprintf(iztex,"\n\\begin{axis}                                                     ");
 fprintf(iztex,"\n[title={ENTER YOUR TITLE HERE},                                   ");  
 fprintf(iztex,"\nwidth=26cm,compat=1.5,height=9cm,                                 ");   
 fprintf(iztex,"\nxtick={1,100,..., 1000},                                          ");
 fprintf(iztex,"\nxticklabels={%i",left_interval);
 for (i=1;i<11;i++)
     {
      fprintf(iztex,",%i",left_interval+(right_interval-left_interval+1)/10*i);
     }
 fprintf(iztex,"},");
 fprintf(iztex,"\nytick=\\empty,yticklabels=\\empty,                                ");       


 fprintf(iztex,"\nxmin=1.0,xmax=1000.0,");
 fprintf(iztex,"\nymin=0,ymax=40,");         
 fprintf(iztex,"\ngrid=both,grid style={line width=.1pt, draw=gray!10},major grid style={line width=.2pt,draw=gray!50},minor tick num=5,enlargelimits={abs=0.5},                  ");
 fprintf(iztex,"\n]                  "); 
 fprintf(iztex,"\n\\draw [green] (1.0,10.0) -- (1000.0,10.0);                  ");    
 
 readLine = fgets(line,101,input);
    while (readLine != NULL && !feof(input))
    {
      sscanf(line,"%[^,],%i,%i,%i,%i,%i,%[^,],%s",protein,&levo1,&desno1,&levo2,&desno2,&duzina,leva_niska,desna_niska);
      if ((left_interval>levo1) || (right_interval<desno2) || ((niska==1) && (strcmp(sequence,leva_niska)!=0)) || ((gel>0) && (duzina<sequence_length)) || ((el>0) && (duzina!=sequence_length))) 
         { 
          /* skip this row */
         }
         else 
          {   
           skalirana_vrednost1=(levo1*1.0-left_interval*1.0+1.0)*1000.0/(right_interval*1.0-left_interval*1.0);
           skalirana_vrednost2=(levo2*1.0-left_interval*1.0+1.0)*1000.0/(right_interval*1.0-left_interval*1.0);
	   fprintf(iztex,"\n\\draw [red,->] (%f,10.0) to [bend left] (%f,10.0); %% row from %s",skalirana_vrednost1, skalirana_vrednost2,line);
          }
      readLine = fgets(line,101,input);
    };
 

fprintf(iztex,"\n\\end{axis}                                ");  
fprintf(iztex,"\n\\end{tikzpicture}                         "); 
fprintf(iztex,"\n\\end{landscape}                           ");
fprintf(iztex,"\n\\end{document}                            ");   
fclose(iztex);

fclose(input);


}