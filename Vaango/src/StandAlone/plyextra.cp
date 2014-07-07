/*

Header for PLY polygon files.

- Greg Turk, March 1994

A PLY file contains a single polygonal _object_.

An object is composed of lists of _elements_.  Typical elements are
vertices, faces, edges and materials.

Each type of element for a given object has one or more _properties_
associated with the element type.  For instance, a vertex element may
have as properties three floating-point values x,y,z and three unsigned
chars for red, green and blue.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

/*
--------------------------------------------------------------------------------
Joao Fradinho Oliveira, July 2005
Copyright (c) 2005 University College London 
copyright conditions as above

update for ply reading of multi OS ply files, in any OS (Unix, Macintosh, PC)
--------------------------------------------------------------------------------

int read_PLY_X(char *fname) 

--------- sets up temp structures
--------- to call ply_open_for_reading library multi OS capable modified function
--------- shows how to release temporary list memory allocated by lib, whilst
--------- inserting read data in external user data structure



int read_PLY_binary(char *fnx, int verts, int ferts)

--------- alternative ply binary read only(not ascii), not multiOS, but fast if you already know what you are reading
--------- copes with little endian/big endian, using functions below



int write_PLY_X(char *fname, int PLY_MODE, int OSlinebreaks=0)

--------- writes ply files either in binary or ascii, with the chosen line breaks convention
--------- PLYMODE 1 (ascii) PLYMODE 2 (binary) 
--------- OSlinebreaks 0 (mac), OSlinebreaks 1 (unix), OSlinebreaks 2 (pc)



--------- binary reading and swapping for big endian and little endian conversion

unsigned char getbyte(FILE *inf)
unsigned int getuint_NATIVE(FILE *inf)
unsigned int getuint_SWAP(FILE *inf)
int getint_NATIVE(FILE *inf)
int getint_SWAP(FILE *inf)
float getfloat_NATIVE(FILE *inf)
float getfloat_SWAP(FILE *inf)
double getdouble_NATIVE(FILE *inf)
double getdouble_SWAP(FILE *inf)
int putbyte(FILE *outf, unsigned char val)
int putshort_NATIVE(FILE *outf, unsigned short val)
int putuint_NATIVE(FILE *outf, unsigned int val)
int putint_NATIVE(FILE *outf, int val)
int putfloat_NATIVE(FILE *outf, float val)
int putdouble_NATIVE(FILE *outf, double val)



*/



/* user's vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;             /* the usual 3-space position of a vertex */
} Vertex;

typedef struct FaceX {
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
} FaceX;


/* information needed to describe the user's data to the PLY routines */

//static char *elem_names[] = { /* list of the kinds of elements in the user's object */
//  "vertex", "face"
//};

static PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,z), 0, 0, 0, 0},
};

static PlyProperty face_props[] = { /* list of property information for a vertex */
  {"vertex_indices", PLY_INT, PLY_INT, offsetof(FaceX,verts),
   1, PLY_UCHAR, PLY_UCHAR, offsetof(FaceX,nverts)},
};




template<class P, class F, class E> int GObject<P,F,E>::read_PLY_X(char *fname)
{
	float fa, fb, fc, fd;
	unsigned int nid;
	
	char str[200];
	int tokencount;
	stringutil temptoks;
	unsigned int countv=0;
	unsigned int countf=0;
	unsigned int intcount=0;
 	unsigned int prevcount;
	unsigned int iii(0);
	unsigned int vstart, fstart;
	GObject<P,F,E> *tempobject=0;

  /////////////////////////////////////
  int i,j,k;
  PlyFile *ply;
  PlyOtherElems *other_elements = NULL;

  int nelems;
  char **elist;
  int file_type;
  float version;
  int nprops;
  int num_elems;
  PlyProperty **plist;
  Vertex **vlist;
  FaceX **flist;
  char *elem_name;
  int num_comments;
  char **comments;
  int num_obj_info;
  char **obj_info;
  /////////////////////////////////////
	
	//this->set_object_name(fname);
	tempobject=new GObject<P,F,E>();
	vstart=this->sizeOfVertexArray();
	fstart=this->sizeOfFaceArray();
	Normals.read_normals(5);
	cout<<"reading..."<<fname<<"\n";

/////////////////////////////////////////////        	
	 /* open a PLY file for reading */
  ply = ply_open_for_reading(fname, &nelems, &elist, &file_type, &version);

  /* print what we found out about the file */
  printf ("version %f\n", version);
  printf ("type %d\n", file_type);

  /* go through each kind of element that we learned is in the file */
  /* and read them */
  
   for (i = 0; i < nelems; i++) {

    // get the description of the first element  
    elem_name = elist[i];
    plist = ply_get_element_description (ply, elem_name, &num_elems, &nprops);

    // print the name of the element, for debugging  
    printf ("element %s %d\n", elem_name, num_elems);

    // if we're on vertex elements, read them in  
    if (equal_strings ("vertex", elem_name)) {

      // create a vertex list to hold all the vertices  
      // vlist = (Vertex **) malloc (sizeof (Vertex *) * num_elems);
      vlist = (Vertex **) malloc (sizeof (Vertex *) * 1); //one vertex only

      // set up for getting vertex elements  

      ply_get_property (ply, elem_name, &vert_props[0]);
      ply_get_property (ply, elem_name, &vert_props[1]);
      ply_get_property (ply, elem_name, &vert_props[2]);

 	  vlist[0] = (Vertex *) malloc (sizeof (Vertex));
      // grab all the vertex elements  
      for (j = 0; j < num_elems; j++) {

	 
        // grab and element from the file  
		
		//vlist[j] = (Vertex *) malloc (sizeof (Vertex));
        ply_get_element (ply, (void *) vlist[0]);

		 this->insertPoint3D(vlist[0]->x, vlist[0]->y, vlist[0]->z);
		//this->insertPoint3D(1,2, 3);

		countv++;

        // print out vertex x,y,z for debugging  
        //printf ("vertex: %g %g %g\n", vlist[0]->x, vlist[0]->y, vlist[0]->z);
      	}
      }
    
      // if we're on face elements, read them in  
    else if (equal_strings ("face", elem_name)) {

      // create a list to hold all the face elements 
      //flist = (Face **) malloc (sizeof (Face *) * num_elems);
      flist = (FaceX **) malloc (sizeof (FaceX *) * 1); // one face


      // set up for getting face elements  

       ply_get_property (ply, elem_name, &face_props[0]);
 
      flist[0] = (FaceX *) malloc (sizeof (FaceX));
      
      Face fjoao;
      
      // grab all the face elements  
      for (j = 0; j < num_elems; j++) {

        // grab and element from the file  
        	(flist[0])->verts=NULL;

         ply_get_element (ply, (void *) flist[0]);
        // print out face info, for debugging  
        //printf ("face: %d, list = ", flist[0]->intensity);
             
        fjoao.set_verts(flist[0]->verts[0], flist[0]->verts[1], flist[0]->verts[2]);
      
          free((flist[0])->verts);	 ////////////// free internal lists allocated by ply library
        
       
       // for (k = 0; k < flist[0]->nverts; k++)
       //   printf ("%d ", flist[0]->verts[k]);
       // printf ("\n");
        
 
         this->insertFace(fjoao);
		countf++; 	
        
      }
    }
    else
    {
		other_elements = ply_get_other_element (ply, elem_name, num_elems);
    }
    
    // print out the properties we got, for debugging  
    for (j = 0; j < nprops; j++)
      printf ("property %s\n", plist[j]->name);
  }    
        
        
        
        
/////////////////////////////////////////////        
	
	 
	tempobject->addFacesToObject(fstart, countf);
	tempobject->addVerticesToObject(vstart, countv);
	this->addChildObject(tempobject);
		
	numVertices=this->sizeOfVertexArray();
	numFaces=this->sizeOfFaceArray();	
		
	
	this->addFacesToObject(0, numFaces);
	this->addVerticesToObject(0, numVertices);

	//tempobject->addfacestovertexneighbors();	


	this->set_PlaneEqOfObject();		
	this->buildnormalsvertex();
	
	cout<<"read "<<countv<<" vertices\n";
	cout<<"read "<<countf<<" faces\n";

 
 
  /////////////////////////////////////////////////////////////////////////////
  /* grab and print out the comments in the file */
  comments = ply_get_comments (ply, &num_comments);
  for (i = 0; i < num_comments; i++)
    printf ("comment = '%s'\n", comments[i]);

  /* grab and print out the object information */
  obj_info = ply_get_obj_info (ply, &num_obj_info);
  for (i = 0; i < num_obj_info; i++)
    printf ("obj_info = '%s'\n", obj_info[i]);

  /* close the PLY file */
  ply_close (ply);
 /////////////////////////////////////////////////////////////////////////////
 
 
 
 	cout<<"reading complete!\n";
	return 1;	
}


//will only read unix line breaks binary files
template<class P, class F, class E> int GObject<P,F,E>::read_PLY_binary(char *fnx, int verts, int ferts)
{
	char str[100];
	int tokencount;
	stringutil temptoks;
	int countv=0;
	int countf=0;
	int intcount=0;
	//time_t time1, time2;
	int prevcount;
	int iii(0);
	int bmode;
	int vstart, fstart;
	GObject<P,F,E> *tempobject=0;

  	fpos_t pos; //keep track of file pointer
  	int nbytes;
  
	FILE *of;
 
 	of = fopen(fnx,"r+b");

	 if(of==NULL) {
	 	fprintf(stderr,"sgiimage: can't open output file\n");
	 	exit(1);
	 }
	
	//	get header, e.g.
	//	ply
	//	format binary_big_endian 1.0
	//  element vertex 3609600
	// 	property float x
	//	property float y
	//	property float z
	//	element face 7219045
	//	property list uchar int vertex_indices
	//	end_header
	
 	fgets(str, 100, of); //depending on line breaks this call can get endian information
	fgets(str, 100, of);
	if(str[14]=='b') //big endian
	{
		bmode=2; //big endian
	}
	else
	{
		bmode=3; //little endian
	}
	fgets(str, 100, of);
	fgets(str, 100, of);
	fgets(str, 100, of);
	fgets(str, 100, of);
	fgets(str, 100, of);
	fgets(str, 100, of);
	 ////fgets(str, 100, of); // get confidence property?
	fgets(str, 100, of);

	
 	tempobject=new GObject<P,F,E>();

	vstart=this->sizeOfVertexArray();
	fstart=this->sizeOfFaceArray();
	
 	int jjj;
 	float vvv1, vvv2, vvv3;
 	unsigned int l;
 	unsigned char uc;
 	
 	F ftemp;
 	double da, db, dc, dd;
 	float fa, fb, fc, fd; 
	Normals.read_normals(5); //read normals
	unsigned short int nid;
 	if(get_native_binary_type2()==bmode) //native binary direct read
	{

		for(jjj=0; jjj<verts; jjj++)
		{

	 		vvv1=getfloat_NATIVE(of);
	 		vvv2=getfloat_NATIVE(of);
	 		vvv3=getfloat_NATIVE(of);

			this->insertPoint3D(vvv1, vvv2, vvv3);
			
			 //vvv3=getfloat_NATIVE(of);//get confidence float property?

			countv++;
 		 
		}
	
	
		for(jjj=0; jjj<ferts; jjj++)
		{
			uc=getbyte(of); //number of vertices is just one byte
		
 			v1idx=getint_NATIVE(of);
 			v2idx=getint_NATIVE(of);
 			v3idx=getint_NATIVE(of);
 		
 			ftemp.set_verts(v1idx, v2idx, v3idx);
 			
			this->insertFace(ftemp);
			
			countf++; 				
		}
  
	}
	else
	{
		for(jjj=0; jjj<verts; jjj++)
		{
 
			
			vvv1=getfloat_SWAP(of);
	 		vvv2=getfloat_SWAP(of);
	 		vvv3=getfloat_SWAP(of);

			this->insertPoint3D(vvv1, vvv2, vvv3);
			
			//vvv3=getfloat_SWAP(of);//get confidence float property?

			countv++;
 		 
		}
	
	
		for(jjj=0; jjj<ferts; jjj++)
		{
			uc=getbyte(of); //number of vertices is just one byte
		
 			v1idx=getint_SWAP(of);
 			v2idx=getint_SWAP(of);
 			v3idx=getint_SWAP(of);
 		
 		
			this->insertFace(v1idx, v2idx, v3idx);
			countf++; 				
		}

	}
 	

	
	tempobject->addFacesToObject(fstart, countf);
	tempobject->addVerticesToObject(vstart, countv);
	this->addChildObject(tempobject);
		
	numVertices=this->sizeOfVertexArray();
	numFaces=this->sizeOfFaceArray();	
		
	
	this->addFacesToObject(0, numFaces);
	this->addVerticesToObject(0, numVertices);

	 tempobject->addfacestovertexneighbors();	

	this->set_PlaneEqOfObject();		
	 this->buildnormalsvertex();
	
	cout<<"read "<<countv<<" vertices\n";
	cout<<"read "<<countf<<" faces\n";

	fclose(of);
	cout<<"reading complete!\n";
	return 1;	
}

template<class P, class F, class E> int GObject<P,F,E>::write_PLY_X(char *fname, int PLY_MODE, int OSlinebreaks=0)
{
	int k, i, j, sz;
	int vcount=0;
	int spec=0;
	
	FILE *out;		// we will use standard C if writing in binary, normal strings are allowed in binary
	ofstream outa;	// we will use a stream if writing in ascii
	
	int v1idx, v2idx, v3idx;
	float x, y, z;
	unsigned int argh;
	unsigned char vn;

	
	P *p;
	F *f;
	
 	// on a macintosh a line break is one character    				 '\r'
	// on a unix machine a line break is one character				 '\n'
	// on a pc a line break is two characters, or the string         "\r\n"
	// note: do not trust programs that convert line breaks to pc, specially if it contains binary
	
	 
	
	//check to see if fname is completely specified, e.g did user provide file extension on the name
	sz=strlen(fname);
	char *fname2= new char [sz+4];
	for(k=0; k<sz; k++)
	{
		if(fname[k]=='.')
		{spec=1;break;}
	}
	strcpy(fname2, fname);
	
	if(spec) //fully specified file
	{;}
	else	//partial filename given, need to add extension
	{
		strcat(fname2, ".ply");
	}
	
	// open in binary, even if writing in ascii
	// in binary mode, there are no translations with line breaks, so we can create
	// line breaks in the OS we choose
	
	out = fopen(fname2,"wb");	
	delete fname2;
	if(!out)
	{
		cout<<"Cannot open file.\n";
		return 1;
	}
		
	// write the header	
	if(OSlinebreaks==0) // Macintosh
	{
		fprintf(out,"ply\r");
		if(PLY_MODE==1)
		{
			fprintf(out,"format ascii 1.0\r");
		}
		else
		{
				if(get_native_binary_type2()==PLY_BINARY_BE)
				{
					fprintf(out,"format binary_big_endian 1.0\r");
				}
				else
				{
					fprintf(out,"format binary_little_endian 1.0\r");
				}
		}	
		fprintf(out,"element vertex %d\r", this->actualSizeOfVertexArray());
		fprintf(out,"property float x\r");
		fprintf(out,"property float y\r");
		fprintf(out,"property float z\r");
		fprintf(out,"element face %d\r", this->actualSizeOfFaceArray());
		fprintf(out,"property list uchar int vertex_indices\r");
		fprintf(out,"end_header\r");
	}
	else if(OSlinebreaks==1) // Unix
	{
		fprintf(out,"ply\n");
		if(PLY_MODE==1)
		{
			fprintf(out,"format ascii 1.0\n");
		}
		else
		{
				if(get_native_binary_type2()==PLY_BINARY_BE)
				{
					fprintf(out,"format binary_big_endian 1.0\n");
				}
				else
				{
					fprintf(out,"format binary_little_endian 1.0\n");
				}
		}
		fprintf(out,"element vertex %d\n", this->actualSizeOfVertexArray());
		fprintf(out,"property float x\n");
		fprintf(out,"property float y\n");
		fprintf(out,"property float z\n");
		fprintf(out,"element face %d\n", this->actualSizeOfFaceArray());
		fprintf(out,"property list uchar int vertex_indices\n");
		fprintf(out,"end_header\n");
	}
	else if(OSlinebreaks==2) // pc
	{
		fprintf(out,"ply\r\n");
		if(PLY_MODE==1)
		{
			fprintf(out,"format ascii 1.0\r\n");
		}
		else
		{
				if(get_native_binary_type2()==PLY_BINARY_BE)
				{
					fprintf(out,"format binary_big_endian 1.0\r\n");
				}
				else
				{
					fprintf(out,"format binary_little_endian 1.0\r\n");
				}
		}
		fprintf(out,"element vertex %d\r\n", this->actualSizeOfVertexArray());
		fprintf(out,"property float x\r\n");
		fprintf(out,"property float y\r\n");
		fprintf(out,"property float z\r\n");
		fprintf(out,"element face %d\r\n", this->actualSizeOfFaceArray());
		fprintf(out,"property list uchar int vertex_indices\r\n");
		fprintf(out,"end_header\r\n");
	}
	
	// write data
	if(PLY_MODE!=1) //not ASCII
	{
			for(j=0; j<this->sizeOfVertexArray(); j++)
			{
				p=this->atVertexArray(j);
				if(p!=0)
				{
					x=p->X(); y=p->Y(); z=p->Z();
					putfloat_NATIVE(out, x);
					putfloat_NATIVE(out, y);
					putfloat_NATIVE(out, z);
 				}
 			}	
	
 			vn=3;
 	
			for(j=0; j<this->sizeOfFaceArray(); j++)
			{
				f=this->atFaceArray(j);
				if(f!=0)
				{
					v1idx=f->get_firstvertexid(); 
					v2idx=f->get_secondvertexid(); 
					v3idx=f->get_thirdvertexid(); 
		 
					putbyte(out, vn);
					putint_NATIVE(out, v1idx); 	
					putint_NATIVE(out, v2idx); 	
					putint_NATIVE(out, v3idx); 	
 				}
			}
	}
	else	//ASCII 1
	{
		if(OSlinebreaks==0) // Macintosh
		{
			for(j=0; j<this->sizeOfVertexArray(); j++)
			{
				p=this->atVertexArray(j);
				if(p!=0)
				{				
						fprintf(out,"%f %f %f\r", p->X(), p->Y(), p->Z());
 				}
 			}
 		
 			for(j=0; j<this->sizeOfFaceArray(); j++)
			{
				f=this->atFaceArray(j);
				if(f!=0)
				{
					v1idx=f->get_firstvertexid(); 
					v2idx=f->get_secondvertexid(); 
					v3idx=f->get_thirdvertexid(); 
		 	
		 			fprintf(out,"3 %d %d %d\r", v1idx, v2idx, v3idx);

				}	 
			}
		}
		else if(OSlinebreaks==1) // Unix
		{
			for(j=0; j<this->sizeOfVertexArray(); j++)
			{
				p=this->atVertexArray(j);
				if(p!=0)
				{				
						fprintf(out,"%f %f %f\n", p->X(), p->Y(), p->Z());
 				}
 			}
 		
 			for(j=0; j<this->sizeOfFaceArray(); j++)
			{
				f=this->atFaceArray(j);
				if(f!=0)
				{
					v1idx=f->get_firstvertexid(); 
					v2idx=f->get_secondvertexid(); 
					v3idx=f->get_thirdvertexid(); 
		 	
		 			fprintf(out,"3 %d %d %d\n", v1idx, v2idx, v3idx);

				}	 
			}
		}
		else if(OSlinebreaks==2)  // pc
		{
			for(j=0; j<this->sizeOfVertexArray(); j++)
			{
				p=this->atVertexArray(j);
				if(p!=0)
				{				
						fprintf(out,"%f %f %f\r\n", p->X(), p->Y(), p->Z());
 				}
 			}
 		
 			for(j=0; j<this->sizeOfFaceArray(); j++)
			{
				f=this->atFaceArray(j);
				if(f!=0)
				{
					v1idx=f->get_firstvertexid(); 
					v2idx=f->get_secondvertexid(); 
					v3idx=f->get_thirdvertexid(); 
		 	
		 			fprintf(out,"3 %d %d %d\r\n", v1idx, v2idx, v3idx);

				}	 
			}
		}			
	
	}
	
	fclose(out); 
		
		
	
 	cout<<"writing complete!\n";
	return 0;
}



template<class P, class F, class E> unsigned char GObject<P,F,E>::getbyte(FILE *inf)
{
	unsigned char buf[1];

	fread(buf, 1, 1, inf);
	
	return(buf[0]);
}

template<class P, class F, class E> unsigned int GObject<P,F,E>::getuint_NATIVE(FILE *inf)
{

 	unsigned int temp;
 	
 	fread(&temp, 4, 1, inf);
	
	return temp;
}


template<class P, class F, class E> unsigned int GObject<P,F,E>::getuint_SWAP(FILE *inf)
{

	unsigned char buf[4]; //read 4 bytes, the first read byte is the high order bits, shift them 16 bit positions
	unsigned char buf2[4]; 
	unsigned int temp;
	fread(&buf, 4, 1, inf);
	buf2[0]=buf[3];
	buf2[1]=buf[2];
	buf2[2]=buf[1];
	buf2[3]=buf[0];

	temp=*((unsigned int *)(&buf2));
	return temp;
}

template<class P, class F, class E> int GObject<P,F,E>::getint_NATIVE(FILE *inf)
{

	// read 4 bytes
 	int temp;

 	fread(&temp, 4, 1, inf);
	
	return temp;
}

template<class P, class F, class E> int GObject<P,F,E>::getint_SWAP(FILE *inf)
{

	unsigned char buf[4]; //read 4 bytes, the first read byte is the high order bits, shift them 16 bit positions
	unsigned char buf2[4]; 
	int temp;
	fread(&buf, 4, 1, inf);
	buf2[0]=buf[3];
	buf2[1]=buf[2];
	buf2[2]=buf[1];
	buf2[3]=buf[0];

	temp=*((int *)(&buf2));
	return temp;

}

template<class P, class F, class E> float GObject<P,F,E>::getfloat_NATIVE(FILE *inf)
{

	// read 4 bytes 
 	float temp;
 
	fread(&temp, 4, 1, inf);
	
	return temp;
}

template<class P, class F, class E> float GObject<P,F,E>::getfloat_SWAP(FILE *inf)
{

	unsigned char buf[4]; //read 4 bytes, the first read byte is the high order bits, shift them 16 bit positions
	unsigned char buf2[4]; 
	float temp;
	fread(&buf, 4, 1, inf);
	buf2[0]=buf[3];
	buf2[1]=buf[2];
	buf2[2]=buf[1];
	buf2[3]=buf[0];

	temp=*((float *)(&buf2));
	return temp;


}

template<class P, class F, class E> double GObject<P,F,E>::getdouble_NATIVE(FILE *inf)
{

	// read 4 bytes 
 	double temp;
 
	fread(&temp, 4, 1, inf);
	
	return temp;
}

template<class P, class F, class E> double GObject<P,F,E>::getdouble_SWAP(FILE *inf)
{

	unsigned char buf[8]; //read 4 bytes, the first read byte is the high order bits, shift them 16 bit positions
	unsigned char buf2[8]; 
	double temp;
	fread(&buf, 8, 1, inf);
	buf2[0]=buf[7];
	buf2[1]=buf[6];
	buf2[2]=buf[5];
	buf2[3]=buf[4];
	buf2[4]=buf[3];
	buf2[5]=buf[2];
	buf2[6]=buf[1];
	buf2[7]=buf[0];

	temp=*((double *)(&buf2));
	return temp;


}

template<class P, class F, class E> int GObject<P,F,E>::putbyte(FILE *outf, unsigned char val)
{
	return fwrite(&val,1,1,outf);
}

template<class P, class F, class E> int GObject<P,F,E>::putshort_NATIVE(FILE *outf, unsigned short val)
{
	return fwrite(&val,2,1,outf);
}

template<class P, class F, class E> int GObject<P,F,E>::putuint_NATIVE(FILE *outf, unsigned int val)
{
	return (fwrite(&val,4,1,outf));

}

template<class P, class F, class E> int GObject<P,F,E>::putint_NATIVE(FILE *outf, int val)
{
	unsigned int val2;
 	val2= (*(unsigned int *) &val);
 	
	return (fwrite(&val2,4,1,outf));

}

template<class P, class F, class E> int GObject<P,F,E>::putfloat_NATIVE(FILE *outf, float val)
{
	unsigned int val2;
 	val2= (*(unsigned int *) &val);
	
	return (fwrite(&val2,4,1,outf));
	
}

template<class P, class F, class E> int GObject<P,F,E>::putdouble_NATIVE(FILE *outf, double val)
{
	unsigned long int val2;
 	val2= (*(unsigned long int *) &val);
	
	return (fwrite(&val2,8,1,outf));
	
}
 




