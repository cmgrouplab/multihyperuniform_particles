//mimic the packing of retena cells using a generalized cherry-pit model
//(1)the inner core of the cherry represents the physical cell boundary, which is impentrable
//(2)the outter shell of the cherry represent soft-repulsive long-range interaction between same specices

//we start with regular triangluar arrays of the outter shells with each species
// first randomly superpose different species
// then (1) remove overlap using MC as a hard sphere system
//      (2) quence the system considering the interaction between the outter shell


//the inner core is then allowed to grow, keeping the size ratio
//after each growth we (1) remove overlap between hard-core using MC
//  (2)quench the soft-core interaction to get an energy minimum



//started: 06/12/1012
//author: Yang Jiao
//email: yjiao@princeton.edu

//modified: 10/24/2023
//allowing specifying soft shelll ranges for individual species

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#define MAXS 5000
#define SD 2 //the space dimension
#define neighbor_max 500 //the maximum number of near neighbor for like-species
double pi = 3.141592654;


int n_s; //the number of different species
int* N_c; //number of cell of each species N_cell[n_s]

double* Rad; //radius of each cell type
double* Range; //the interaction range of like-species
double Rad_max; //the maximum radius among all species, for grid list set-up

int N_tot; //the total number of cells
int* Type_c; //type of each cell Type_cell[N_tot]

double growth_rate; //the rescale rate of cell radius R = (1+gr)*R
double rescale = 0.75; //if the growth cause overlap, rescale the growth rate by a half

double box_length = 1.0; //the length of the simulation box

double** coords; //the coords of cell centers

double neighbor_cut = 1.5; //1.5*Range for the like-species
int** NNL; //nnl for the like-species, updated in the beginning of each quench
//allow randomness in the configuration
int* NNL_counter; //the number of nn of like species for each cell

double old_coord[SD]; //the old coord of a moved cell
int move_index; //the index of the moved cell

double trans_mod;

int N_stage; //the total number of particle growth stages...

int N_relax; //the total number of energy minimization relax

double rho_start; //the packing density at which start quenching....

int flag_config = 0; //to generate new initial config (0) or to read old ones (1)

int N_MC_particle = 200; //number of MC moves per particle

double d0_cell = 0.12; //the length of the grid cell



//the following is for the grid (cell) method
//##################################################

double rho_cell; //the density of last update point

//the number of cells along each direction, determined automatically later
//    these are also updated adaptively
int LXcell = 3; 
int LYcell = 3;
int Ncell = LXcell*LYcell; //the total number of cells

struct node
{
  int index_poly;
  node* next;
};

node** CellList; //the cell list containing the particles in the cell...
//the corner locations are the index...
//implemented as a one-dimensional array




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//functions related to the cell method

int NumMax(int a, int b)
{
  if(a>b) return a;
  else return b;
}

void DelCellList()
{
  for(int i=0; i<Ncell; i++)
    { 
      while(CellList[i]!=NULL)
	{
	  node* temp_pt = CellList[i];
	  CellList[i] = CellList[i]->next;
	  
	  delete temp_pt;
	}
    }
}


void GetCellList(int flag, double &d0, double** &CenterL, int &N, double &temp_rho)
{
  
  if(flag != 0) //not the first time to establish CellList, so need to clean the list first
    {
      DelCellList();
    }

  //now need to compute the number of cell along each direction
  double LatX, LatY; //the length of lattice vectors along each direction
  
  LatX = box_length;
  LatY = box_length;

  //double linear_size = sqrt((LatX*LatY)/(double)N);

  LXcell = (int)NumMax(3, (int)floor(LatX/(1.25*d0)));
  LYcell = (int)NumMax(3, (int)floor(LatY/(1.25*d0)));
  //the size of cell has to be large enough such that next neighbor cell can never contain possibly overlapped particles

  Ncell = LXcell*LYcell;
  
  cout<<"The total number of cells Ncell = "<<Ncell<<endl;

  CellList = new node*[Ncell]; //re-setup the cell list
  for(int i=0; i<Ncell; i++)
    CellList[i] = NULL;

  for(int n=0; n<N; n++)
    {
      int temp_indexI = (int)floor(CenterL[n][0]*LXcell);
      int temp_indexJ = (int)floor(CenterL[n][1]*LYcell);

      if(temp_indexI>=LXcell || temp_indexJ>=LYcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  cout<<"cell index n = "<<n<<endl;
	  cout<<CenterL[n][0]<<" "<<CenterL[n][1]<<endl;
	  exit(1);
	}

      int cell_index = LXcell*temp_indexJ + temp_indexI;

      //cout<<"Tetrah"<<n<<" cell_index = "<<cell_index<<endl;

      node* pt = new node;
      pt->index_poly = n;
      pt->next = CellList[cell_index];
      CellList[cell_index] = pt;

    }

  //record the density of the most recent update
  rho_cell = temp_rho;
  
}


void UpdateCellList(int m, double Old_CenterL[SD], double CenterLm[SD])
{
  //delete the particle from the old list...
  int temp_indI = (int)floor(Old_CenterL[0]*LXcell);
  int temp_indJ = (int)floor(Old_CenterL[1]*LYcell);

  int cell_index = LXcell*temp_indJ + temp_indI;

  node* pt1 = CellList[cell_index];
  node* pt2 = CellList[cell_index];

  if(pt1 == NULL)
    { 
       printf("Bugs exist in CellList! Re-check!\n");
       //print out the list..
       printf("The particle is %d\n", m);
       printf("It should be in (%d, %d)\n", temp_indI, temp_indJ);
       //PrintCellList();
    
       exit(1);
    }
  else
    {
      while(pt1!=NULL)
	{
	  if((pt1->index_poly)==m)
	    {
	      if(pt2==CellList[cell_index] && pt1==CellList[cell_index])
		CellList[cell_index] = pt1->next;
	      else
                pt2->next = pt1->next;
	      free(pt1);
	      break;
	    }

	  pt2 = pt1;
	  pt1 = pt1->next;

	}
      /*
      if(pt1==NULL)
	{
	  printf("Bugs exist in celllist! Recheck!\n");
	  exit(1);
	}
      */
    }


  //insert the particle in the new list...
  temp_indI = (int)floor(CenterLm[0]*LXcell);
  temp_indJ = (int)floor(CenterLm[1]*LYcell);

  if(temp_indI>=LXcell || temp_indJ>=LYcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  cout<<"cell index n = "<<m<<endl;
	  cout<<CenterLm[0]<<" "<<CenterLm[1]<<endl;
	  exit(1);
	}

  int new_cell_index = LXcell*temp_indJ + temp_indI;

  node* pt = (node *)malloc(sizeof(node));
  pt->index_poly = m;
  pt->next = CellList[new_cell_index];

  CellList[new_cell_index] = pt;

}

void PrintCellList()
{
   for(int i=0; i<Ncell; i++)
     {
       node* pt = CellList[i];

       int temp_indexJ = (int)floor((double)i/(double)LXcell);
       int temp_indexI = (i%LXcell);
      
       printf("CellList(%d,%d) = ", temp_indexI, temp_indexJ);
       
       while(pt!=NULL)
	 {
	   printf(" %d ", pt->index_poly);
	   pt = pt->next;
	 }
       
       printf("\n");
       
     }
   
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void read_parameter()
{
  cout<<"total number of different species n_s = "; cin>>n_s;
  N_c = new int[n_s];
  Rad = new double[n_s];
  Range = new double[n_s];
  N_tot = 0;
  N_stage = 0;
  N_relax = 0;

  //this is only a rough number determining the long-interaction range
  //tha actual number is detemined by how many can fit into the unit box on a triangular lattice
  
  cout<<"the total number of growth stages N_stage = "; cin>>N_stage; cout<<endl;

  cout<<"the total number of relax stages N_relax = "; cin>>N_relax; cout<<endl;

  for(int i=0; i<n_s; i++)
    {
      cout<<"number of cells for species "<<i<<" N_c["<<i<<"] = "; cin>>N_c[i];
      N_tot += N_c[i];
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  for(int i=0; i<n_s; i++)
    {
      cout<<"radius of cells for species "<<i<<" Rad["<<i<<"] = "; cin>>Rad[i];
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  for(int i=0; i<n_s; i++)
    {
      Range[i] = sqrt(box_length*box_length/(pi*N_c[i]));
      //this is the RSA density 
      cout<<"MAXIMAL interfaction range of cells for species "<<i<<" Range["<<i<<"] = "<<Range[i]<<endl;
      cout<<"Rescale factor (preferred to be between 0 and 2, larger values might cause error) for species" <<i<<endl;
      
      //now allow rescaling the soft shell
      double tmp_range_scale;
      cin>>tmp_range_scale;
      
      Range[i] = tmp_range_scale*Range[i];
      cout<<"RESCALED interfaction range of cells for species "<<i<<" Range["<<i<<"] = "<<Range[i]<<endl;
      
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"the growth rate (1+gr)*Rad growth_rate = "; cin>>growth_rate;

  cout<<"the trans_mod (in terms of boxlength) ="; cin>>trans_mod;

  cout<<"the start quenching density rho_start = "; cin>>rho_start;

  cout<<"generate initial config or read in flag_config = "; cin>>flag_config; 

  cout<<"number of MC moves per particle at each stage N_MC_particle = "; cin>>N_MC_particle;

  cout<<"the length of the grid cell: d0_cell = "; cin>>d0_cell;

  //now we initialize the quantities with the cell information
  N_tot = 2*N_tot; //this is a safe guard, when generating configuration, the true N_tot would be larger than this value

  coords = new double*[N_tot];
  for(int i=0; i<N_tot; i++)
    {
      coords[i] = new double[SD];
    }

  NNL = new int*[N_tot];
  for(int i=0; i<N_tot; i++)
    {
      NNL[i] = new int[neighbor_max];
    }

  NNL_counter = new int[N_tot];

  Type_c = new int[N_tot];
  


}


double dist(int index1, int index2)
{
  double tempx = fabs(coords[index1][0] - coords[index2][0]);
  if(tempx >= box_length/2.0) tempx = box_length - tempx;

  double tempy = fabs(coords[index1][1] - coords[index2][1]);
  if(tempy >= box_length/2.0) tempy = box_length - tempy;

  return sqrt(tempx*tempx+tempy*tempy);

  
}


double dist(double center1[SD], double center2[SD])
{
  double tempx = fabs(center1[0] - center2[0]);
  if(tempx >= box_length/2.0) tempx = box_length - tempx;

  double tempy = fabs(center1[1] - center2[1]);
  if(tempy >= box_length/2.0) tempy = box_length - tempy;

  return sqrt(tempx*tempx+tempy*tempy);

  
}


void read_config()
{
  N_tot = 0;
  int temp_Nc = 0;
  int counter = 0;
  int temp_N_tot;
  int counter_global = 0;

  ifstream fin;
  fin.open("cell_config.txt");
  
  fin>>temp_N_tot;

  for(int s=0; s<n_s; s++)
    {
      fin>>temp_Nc;
      cout<<"number of "<<s<<" = "<<temp_Nc<<endl;
      N_c[s] = temp_Nc;

      fin>>Rad[s]; cout<<"Rad["<<s<<"] = "<<Rad[s]<<endl;
      fin>>Range[s];

      for(int i=0; i<N_c[s]; i++)
	{
	  fin>>coords[counter][0];
	  fin>>coords[counter][1];

	  counter++;


	  //now specify the cell types
	  Type_c[counter_global] = s;
	  counter_global++;
	  
	}

      N_tot += N_c[s];
    
      
    }

  fin.close();
  
  if(temp_N_tot != N_tot)
    {
      cout<<"N_tot is not set correctly"<<endl;
      exit(1);
    }
  
}


//use a 2-scale RSA to generate the initial configuration

void get_config()
{
  
  int counter = 0; //this is for each species
  int start = 0; //this is for the outter-shell overlapping check
  N_tot = 0;
  
  double temp_center[2];


  cout<<"generating initial configurations..."<<endl;


  for(int s=0; s<n_s; s++)
    {
      counter = 0;
      start = 0;
      
      for(int ts = 0; ts<s; ts++)
	{
	  start += N_c[ts];
	}

      
      for(int i=0; i<N_c[s]; i++)
	{
	  
	  int overlap_flag = 0;
	  int limit = 0;

	  temp_center[0] = (double)(rand()%MAXS)*box_length/(double)MAXS;
	  temp_center[1] = (double)(rand()%MAXS)*box_length/(double)MAXS;

	  double temp_dist;

	  //start looping over all particles
	  for(int j=0; j<N_tot; j++)
	    {
	      temp_dist = dist(temp_center, coords[j]);

	      //first, check the un-like species
	      if(temp_dist<(Rad[Type_c[j]]+Rad[s]))
		{
		  overlap_flag = 1;
		  break;
		}
	      else
		{
		  //now, check the like species

		  if(j >= start)
		    {
		      	
		      if(temp_dist<0.5*Range[s] || temp_dist<2*Rad[s]) //this is to make sure the like-species do not get too close, but this condition is too strong for RSA!
			  {
			    overlap_flag = 1;
			    break;
			  }
		      
		      
		      /*
		      if(temp_dist<2*Rad[s]) //this is to make sure the like-species do not overlap by their hard core, may cause high energy for larage soft shell radius
			  {
			    overlap_flag = 1;
			    break;
			  }
		      */
		    }
		}
	    }
	  //end looping over all particles


	  while(overlap_flag == 1 && limit < 3000)
	    {
	      overlap_flag = 0;

	      temp_center[0] = (double)(rand()%MAXS)*box_length/(double)MAXS;
	      temp_center[1] = (double)(rand()%MAXS)*box_length/(double)MAXS;

	      double temp_dist;

	      //start looping over all particles
	      for(int j=0; j<N_tot; j++)
		{
		  temp_dist = dist(temp_center, coords[j]);
		  
		  //first, check the un-like species
		  if(temp_dist<(Rad[Type_c[j]]+Rad[s]))
		    {
		      overlap_flag = 1;
		      break;
		    }
		  else
		    {
		      //now, check the like species
		      
		      if(j >= start)
			{
				
			  if(temp_dist<0.5*Range[s] || temp_dist<2*Rad[s]) //this is to make sure the like-species do not get too close, but this condition is too strong for RSA!
			    {
			      overlap_flag = 1;
			      break;
			    }
			    
			    
			    /*
			      if(temp_dist<2*Rad[s]) //this is to make sure the like-species do not overlap by their hard core, may cause high energy for larage soft shell radius
			   {
			    overlap_flag = 1;
			    break;
			   }
		      */
			  
			}
		    }
		}
	     

	      limit++;
	    }
	  

	  if(limit<3000)
	    {
	      coords[N_tot][0] = temp_center[0];
	      coords[N_tot][1] = temp_center[1];
	      
	      Type_c[N_tot] = s;
	      
	      counter++;
	      N_tot++;
	      
	      //cout<<"N_tot = "<<N_tot<<endl;
	      
	    }
	  //otherwise, do not insert a new cell of this species


	}
      
       N_c[s] = counter;
       cout<<"N_c["<<s<<"] = "<<counter<<endl;

    }

  
  //after the RSA, we re-set the radius of the outter shell (Range) so that the current configuration 
  // corresponds to an overlapping configuration of the outter shell...
  //this is dis-abled for the control-range simulations....
  /*
  for(int i=0; i<n_s; i++)
    {
      Range[i] = 3.5*sqrt(box_length*box_length/(sqrt(12.0)*N_c[i]));
      //this is the RSA density 
      cout<<"Re-setted interfaction range of cells for species "<<i<<" Range["<<i<<"] = "<<Range[i]<<endl;
    }
  */

}


//will modify this to allow printing out intermediate configurations
void print_config()
{
  ofstream fout;
  fout.open("centers.xls");
  
  int counter = 0;
  int sum = 0;

  fout<<N_tot<<endl<<endl;

  for(int s=0; s<n_s; s++)
    {
      counter = 0;

      fout<<N_c[s]<<endl;
      fout<<Rad[s]<<endl;
      fout<<Range[s]<<endl;
      
      for(int i=sum; i<N_tot; i++)
	{
	  if(Type_c[i] == s)
	    {
	      fout<<coords[i][0]<<"\t"<<coords[i][1]<<endl;

	      counter++;
	    }

	  if(counter == N_c[s])
	    {
	      sum += counter;
	      break;
	    }
	}

      fout<<"//###########################"<<endl;
      fout<<endl<<endl;

    }

  fout.close();
}

void print_hard_config()
{
  ofstream fout;
  fout.open("centers_hard.xls");
  
  int counter = 0;
  int sum = 0;

  fout<<N_tot<<endl<<endl;

  for(int s=0; s<n_s; s++)
    {
      counter = 0;

      fout<<N_c[s]<<endl;
      fout<<Rad[s]<<endl;
      fout<<Range[s]<<endl;
      
      
      for(int i=sum; i<N_tot; i++)
	{
	  if(Type_c[i] == s)
	    {
	      fout<<coords[i][0]<<"\t"<<coords[i][1]<<endl;

	      counter++;
	    }

	  if(counter == N_c[s])
	    {
	      sum += counter;
	      break;
	    }
	}

      fout<<"//###########################"<<endl;
      fout<<endl<<endl;

    }

  fout.close();
}



//update after each radius growth
void get_NNL()
{
  int counter = 0;
  int sum = 0;
  int nnl_ct =0 ;

  double temp_coord[N_tot][SD];
  int temp_index[N_tot];

  for(int s=0; s<n_s; s++)
    {
      counter = 0;
      
      for(int i=sum; i<N_tot; i++)
	{
	  if(Type_c[i] == s)
	    {
	      temp_coord[counter][0] = coords[i][0];
	      temp_coord[counter][1] = coords[i][1];
	      	      
	      temp_index[counter] = i;

	      counter++;
	    }

	  if(counter == N_c[s])
	    {
	      sum += counter;
	      break;
	    }
	}

      for(int i=0; i<counter; i++)
	{
	  nnl_ct = 0;

	  for(int j=0; j<counter; j++)
	    {
	      if(dist(temp_index[i], temp_index[j])<= neighbor_cut*2*Range[s] && temp_index[i]!=temp_index[j])
		{
		  NNL[temp_index[i]][nnl_ct] = temp_index[j];
		  nnl_ct++;

		  if(nnl_ct>=neighbor_max)
		    {
		      cout<<"the number of neighbor is great that neighbor_max!"<<endl;
		      exit(1);
		    }
		}
	    }
	  
	  NNL_counter[temp_index[i]] = nnl_ct;
	}
      
    }

}

//this is for setting up the cell list
void get_Rad_max()
{
  double temp_max = 0.0;

  for(int i=0; i<n_s; i++)
    {
      if(temp_max <= Rad[i])
	temp_max = Rad[i];
    }

  Rad_max = temp_max;
}



int check_overlap(int index, double center[SD])
{
  //delete the particle from the old list...
  int temp_indI = (int)floor(center[0]*LXcell);
  int temp_indJ = (int)floor(center[1]*LYcell);
  
  int cell_index;
  node* pt1;

  for(int i=-1; i<=1; i++)
    for(int j=-1; j<=1; j++)
      {
	int t_ind1 = temp_indI + i;
	if(t_ind1<0) t_ind1 = t_ind1 + LXcell;
	else if(t_ind1>=LXcell) t_ind1 = t_ind1 - LXcell;

	int t_ind2 = temp_indJ + j;
	if(t_ind2<0) t_ind2 = t_ind2 + LYcell;
	else if(t_ind2>=LYcell) t_ind2 = t_ind2 - LYcell;

	cell_index = LXcell*t_ind2 + t_ind1;
	pt1 = CellList[cell_index];

	while(pt1 != NULL)
	  {
	    int cell_index = pt1->index_poly;

	    if(cell_index != index)
	      {
		if(dist(cell_index, index)<=(Rad[Type_c[cell_index]]+Rad[Type_c[index]]))
		  return 1;
	      }

	    pt1 = pt1->next;
	  }
      }

  return 0;
 
 
}

int global_overlap_check()
{
  int overlap_flag = 0;

  for(int i=0; i<N_tot; i++)
    {
      overlap_flag = check_overlap(i, coords[i]);

      if(overlap_flag == 1)
	return overlap_flag;
    }

  return overlap_flag;
  
}

//if move is accepted, may need to update cell_list
int move(int index)
{
  double dx = trans_mod*box_length*((rand()%MAXS/(double)MAXS)-0.5);
  double dy = trans_mod*box_length*((rand()%MAXS/(double)MAXS)-0.5);

  old_coord[0] = coords[index][0];
  old_coord[1] = coords[index][1];
  
  coords[index][0] = coords[index][0] + dx;
  if(coords[index][0]<0) coords[index][0] = coords[index][0] + box_length;
  else if(coords[index][0]>=box_length) coords[index][0] = coords[index][0] - box_length;

  
  coords[index][1] = coords[index][1] + dy;
  if(coords[index][1]<0) coords[index][1] = coords[index][1] + box_length;
  else if(coords[index][1]>=box_length) coords[index][1] = coords[index][1] - box_length;

  // cout<<dx<<", "<<dy<<endl;
  //cout<<coords[index][0]<<", "<<coords[index][1]<<endl;

  move_index = index;

  if(check_overlap(index, coords[index])==1)
    {
      coords[index][0] = old_coord[0];
      coords[index][1] = old_coord[1];

      return 0;
    }

  int old_indI = (int)floor(old_coord[0]*LXcell);
  int old_indJ = (int)floor(old_coord[1]*LYcell);

  int indI = (int)floor(coords[index][0]*LXcell);
  int indJ = (int)floor(coords[index][1]*LYcell);

  if(old_indI != indI || old_indJ != indJ)
    {
      UpdateCellList(index, old_coord, coords[index]);
      
    }
  
  return 1;
}

//this is for the long-range acceptance rule
void resume_move()
{
  int old_indI = (int)floor(old_coord[0]*LXcell);
  int old_indJ = (int)floor(old_coord[1]*LYcell);
  
  int indI = (int)floor(coords[move_index][0]*LXcell);
  int indJ = (int)floor(coords[move_index][1]*LYcell);
  
  if(old_indI != indI || old_indJ != indJ)
    {
      UpdateCellList(move_index, coords[move_index], old_coord);
    }

  coords[move_index][0] = old_coord[0];
  coords[move_index][1] = old_coord[1];
}

int remove_overlap()
{
  //randomly move each particls unitl a overlap-free configuration is obtained
  int limit = 0;
  int overlap_flag = global_overlap_check();
  int move_flag;
  
  while(overlap_flag == 1&& limit<100)
    {
      
      for(int i=0; i<N_tot; i++)
	{
	  int index = i;
	  //cout<<temp_index<<endl;


	  int temp_indI = (int)floor(coords[i][0]*LXcell);
	  int temp_indJ = (int)floor(coords[i][1]*LYcell);
  
	  int cell_index;
	  node* pt1;

	  for(int i=-1; i<=1; i++)
	    for(int j=-1; j<=1; j++)
	      {
		int t_ind1 = temp_indI + i;
		if(t_ind1<0) t_ind1 = t_ind1 + LXcell;
		else if(t_ind1>=LXcell) t_ind1 = t_ind1 - LXcell;
		
		int t_ind2 = temp_indJ + j;
		if(t_ind2<0) t_ind2 = t_ind2 + LYcell;
		else if(t_ind2>=LYcell) t_ind2 = t_ind2 - LYcell;
		
		cell_index = LXcell*t_ind2 + t_ind1;
		pt1 = CellList[cell_index];
		
		while(pt1 != NULL)
		  {
		    int cell_index = pt1->index_poly;
		    
		    if(cell_index != index)
		      {
			if(dist(cell_index, index)<=(Rad[Type_c[cell_index]]+Rad[Type_c[index]]))
			  {
			    cout<<cell_index<<", "<<index<<endl;
			    cout<<dist(cell_index, index)<<"<"<<(Rad[Type_c[cell_index]]+Rad[Type_c[index]])<<endl;

			    
			    move_flag = 0;
			    int temp_limit = 0;
			    while(move_flag == 0 && temp_limit<100)
			      { 
				move(index);
				temp_limit++;
				//cout<<move_flag;
			      }
			    
			  }
		      }
		    
		    pt1 = pt1->next;
		  }
	      }
	  
	  
	}
      cout<<endl;
      


      overlap_flag = global_overlap_check();

      limit++;
      //cout<<"limit = "<<limit<<endl;
    }

  if(limit == 1000)
    {
      //cout<<"can not remove overlap from the current config!"<<endl;
      return 0;
    }
  else
    return 1;
}


double get_density()
{
  double rho = 0.0;

  for(int s=0; s<n_s; s++)
    {
      rho += 3.1415926*Rad[s]*Rad[s]*N_c[s];
    }

  return rho;
}


/*
double get_energy(int index, double center[SD])
{
  double ener = 0;

  for(int i=0; i<NNL_counter[index]; i++)
    {
      if(NNL[index][i] != index)
	{
	  double temp_dist = dist(center, coords[NNL[index][i]]) - 2*Range[Type_c[index]];

	  if(temp_dist<0)
	    {
	      ener += temp_dist*temp_dist;
	    }
	}
    }

  return ener;
}
*/


//loop over all particles of the same type
double get_energy(int index, double center[SD])
{
	double ener = 0;
	
	for(int i=0; i<N_tot; i++)
    {
      if(i != index && Type_c[i] == Type_c[index])
	{
	  double temp_dist = dist(center, coords[i]) - 2*Range[Type_c[index]];

	  if(temp_dist<0)
	    {
	      ener += temp_dist*temp_dist;
	    }
	}
    }

  return ener;
}


double get_cost(int index, double center[SD], double center_old[SD])
{
  double cost = get_energy(index, center) - get_energy(index, center_old);

  return cost;
}

double get_total_ener()
{
	double en_tot = 0;
	for(int i=0; i<N_tot; i++)
    en_tot +=  get_energy(i, coords[i]);
    
    return en_tot;
}

main()
{

  
  srand(time(NULL));
  cout<<"test"<<endl;

  read_parameter();

  if(flag_config == 0)
    get_config();
  else
    read_config();


  //print_config();

  get_NNL();

  get_Rad_max();

  double temp_rho = get_density();

  GetCellList(0, d0_cell, coords, N_tot,  temp_rho);

  PrintCellList();

  remove_overlap();

  print_config();


  double en_int = 0;
  for(int i=0; i<N_tot; i++)
    en_int +=  get_energy(i, coords[i]);
  cout<<"the initial energy of the system E_int = "<<en_int<<endl;
  
  //open up a file to save the density rho at each stage..
  ofstream fout1;
  fout1.open("rho.xls"); //this is to open an empty file
   


  //************************************
  //now start growing the particles and queching...

  int overlap_flag = 0;
  double rho_i = 0;

  int flag_print = 0;
  int relax_counter = 0;

  for(int i=0; i<N_stage; i++)
    {
      //first grow the core and remove overlaps

      for(int s=0; s<n_s; s++)
	{
	  Rad[s] = (1+growth_rate)*Rad[s];
	}
      
      overlap_flag = global_overlap_check();

      int limit = 0;

      double temp_growth_rate = growth_rate;
      
      while(overlap_flag == 1 && limit<50)
	{
	  for(int s=0; s<n_s; s++)
	    {
	      Rad[s] = Rad[s]/(1+temp_growth_rate);
	    }

	  temp_growth_rate = temp_growth_rate*rescale;

	  for(int s=0; s<n_s; s++)
	    {
	      Rad[s] = Rad[s]*(1+temp_growth_rate);
	    }
	  
	  overlap_flag = global_overlap_check();
	}

      if(limit == 50)
	{
	  cout<<"unable to grow the particle at stage  "<<i<<endl;

	  //print_config();

	  //exit(1);
	}
      else
	{
	  //if obtain a configuration with higher density, 
	  //equilibrate the hard-particle configuration first
	  int N_acc = 0;
	  double p_acc = 0.0;
	  int temp_index;

	  int N_MC_hard = N_MC_particle; //number of MC moves per particle
	  int N_MC_soft = N_MC_particle;
	  
	  for(int j=0; j<N_tot*N_MC_hard; j++)
	    {
	      temp_index = rand()%N_tot;

	      N_acc += move(temp_index);
	      	      
	    }

	  p_acc = (double)N_acc/(double)(N_tot*N_MC_hard);
	  cout<<"hard-particle MC p_acc = "<<p_acc<<endl;
	  

	  /*
	  if(p_acc<0.1)
	    {
	      trans_mod = rescale*trans_mod;
	    }
	  */
	  

	  rho_i = get_density();
	  
	 


	  //then quech the configuration to mimize the outter core overlap...
	  if(rho_i > rho_start)
	    {

	      if(relax_counter>N_relax)
		{
		  goto L1;
		}

	      if(flag_print == 0)
		{
		  print_hard_config(); 

		  flag_print = 1;
		  //print the configuration of the hard particle packing after equilibration
		}

		
		cout<<"start relaxing the configuration to minimize E ..."<<endl;

	      relax_counter++;
	      N_acc = 0;
	      double d_E, E_o, E_n;
	      int move_flag;
	      
	      for(int j=0; j<N_tot*N_MC_soft; j++)
		{
		  temp_index = rand()%N_tot;
		  
		  //E_o = get_total_ener();
		  
		  move_flag = move(temp_index);
		  //cout<<"move_flag = "<<move_flag<<endl;
		  
		  while(move_flag == 0) 
		   {
		   	 //resume_move();
		   	 
		   	 //temp_index = rand()%N_tot;
		   	 
		      move_flag = move(temp_index); //make sure the move is success
		   }
		  
		  //E_n = get_total_ener();
		  
		  //now compute the energy cost...
		  
		  d_E = get_cost(temp_index, coords[temp_index], old_coord);
		  //d_E = E_n - E_o;
		  //cout<<"dE ="<<d_E<<endl;
		  
		  if(d_E>0)
		    {
		      resume_move();
		      
		      //en_int = E_o;
		    }
		  else
		    {
		      N_acc++;
		      
		      //en_int = E_n;
		      
		      /*
		      en_int = en_int + d_E;
		      
		      if(en_int < 0) en_int = 0; //overlapping free configuration for soft shell
		      */
		    }
		    
		    //output the total energy for each relaxation step, for minitoring and debugging
		     //cout<<"Relax Step "<<j<<" total energy E ="<<get_total_ener()<<endl;
		  
		}
	      
	      p_acc = (double)N_acc/(double)(N_tot*N_MC_soft);
	      cout<<"sof-core queching p_acc = "<<p_acc<<endl;  
	      
	      //cout<<"final energy of stage "<<i<<" E = "<<en_int<<endl;
	      cout<<"total energy E ="<<get_total_ener()<<endl;

	      print_config();
	    }
	  
	  cout<<"Rad[0] = "<<Rad[0]<<endl;
	  cout<<"rho_["<<i<<" = "<<rho_i<<endl;
	  
	  //now save the density of stage i to file
	  fout1<<i<<"\t"<<rho_i<<endl;
	  
	  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	  
	}
      
      if(i%20 ==0)
      {
	  
       cout<<"Printing the intermediate configuration now for stage "<<i<<endl;
       print_config();
      }

    }


 L1:

  for(int s=0; s<n_s; s++)
    {
      cout<<"rad["<<s<<"] = "<<Rad[s]<<endl;
    }
      
  cout<<"packing density rho = "<<get_density()<<endl;

  //the final configuration is printed
  cout<<"Printing the final configuration now...."<<endl;
  print_config();
  
  fout1.close(); 


  
}
