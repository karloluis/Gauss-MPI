
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "mpi.h"
#include <omp.h>
using namespace std;

/*****************************FUNCTION PROTOTYPES******************************/

int inverse(int, int);
void find_inverses(int, int[]);
void swap(int*, int*, int);


class Space{
  public:
	int p;					//Prime defining Zp
	int* mult_inv;			//Array of multiplicative inverses in Zp
	int* add_inv;			//Array of additive inverses in Zp

	Space(int);
};

class Matrix {
  public:
	int **aug_matrix;		//Matrix being reduced
	int size;				//Number of rows

	void init_with_file(string, int);
	void display();
	void gauss(Space);
	void Resolve(Space);

};

void eliminate(int*, int*, int, int, Space);


/************************************MAIN**************************************/

int main (int argc, char **argv ) {
	int *pivot_row, *row;
	string filename;

	int numprocs, myid;

	MPI::Init(argc,argv);
	numprocs = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();


	int p = 7;				//Prime used to define Zp
	Space modp(p);

		Matrix system;
		system.init_with_file("coefficients.txt", p);

		if(myid == 0){
			cout << endl << system.size << endl << endl;
			system.display();
		}

		//system.gauss(modp);				//implementacion anterior
		system.Resolve(modp);    //aun no acabado

		if(myid == 0){
		system.display();
		}

	MPI::Finalize();

	return 0;
}

/****************************FUNCTION DEFINITIONS******************************/

int inverse(int p, int num) {
//Receives a prime p and num < p, finds the multiplicative inverse of num in Zp
	if (num == 1)
		return 1;

	if (p%num == 0) {
		cerr << "\nError: \n\tProvided value is not a prime number.\n" << endl;
		exit(1);
	}

	else {
		int n = num - inverse(num, p%num);
		return (n*p + 1)/num;
	}
}

void find_inverses(int p, int inv_arr[]) {
	#pragma omp parallel for
	for (int i = 1; i < p; i++)
		inv_arr[i] = inverse(p, i);
	return;
}

void swap(int* row_a, int* row_b, int size){
	int temp;
	#pragma omp parallel for
	for (int i = 0; i < size; i++){
		temp = row_a[i];
		row_a[i] = row_b[i];
		row_b[i] = temp;
	}
}

//Recivimos el pivot row desde broadcast y el row su size y el modp por Send
void eliminate(int* pivotrow, int *row, int pivot, int size, Space modp) {
	int n = modp.add_inv[row[pivot]%modp.p];
	#pragma omp parallel for
	for (int i = 0; i < size + 1; i++) {
		row[i] = row[i] + (pivotrow[i] * n) ;
		row[i] = row[i] %modp.p;
	}
}

Space::Space(int n) {
		p = n;
		int inv;
		//Resize arrays to size p
		mult_inv = new int[p];
		add_inv = new int[p];

		//Fill in the array of multiplicative inverses
		find_inverses(n, mult_inv);
		//Fill in the array of additive inverses
		add_inv[0] = 0;
		//#pragma omp parallel for
		for (int i = 1; i < p/2 + 1; i++) {
			inv = p - i;
			add_inv[i] = inv;
			add_inv[inv] = i;
		}
}

void Matrix::display() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size + 1; j++)
			cout << aug_matrix[i][j] << "\t";
		cout << endl;
	}
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
	return;
}

void Matrix::init_with_file(string filename, int prime) {

	ifstream infile;
	int temp_val, rows = 0, take;
	string temp_line;
	bool tokenerror;

	infile.open(filename.c_str(), fstream::in);
	if (infile.fail()) {
		cerr << "\nError: \n\tFile could not be opened\n" << endl;
		exit(1);
	}

	while (getline (infile, temp_line)) {
		rows++;
	}

	infile.close();

	infile.open(filename.c_str(), fstream::in);
	if (infile.fail()) {
		cerr << "\nError: \n\tFile could not be opened\n" << endl;
		exit(1);
	}

	//Allocate memory and fill in the array
	aug_matrix = new int*[rows];
	for (int i = 0; i < rows; i++){
		aug_matrix[i] = new int[rows + 1];
		for(int j = 0; j < rows + 1; j++){
			infile >> temp_val;
			aug_matrix[i][j] = temp_val % prime;
		}
	}

	size = rows;
	infile.close();
	return;
}



void Matrix::Resolve(Space modp){

	int numprocs, myid, tag = 0;

    numprocs = MPI::COMM_WORLD.Get_size();
    myid = MPI::COMM_WORLD.Get_rank();

	int temp_array[size + 1], *temp_array_pivot = new int[size + 1];
	int i, j, k, inv, j_row;

	int size_chunck, size_chuck_tail ;
	//int temp_chunck*[size+1];

	    if (myid != numprocs -1){
			//chunksize = size/numprocs; 				// chunksize for all node minus the last
	    }
		else if(myid==numprocs-1){
			//chunksize = range - (myid*size/numprocs)+2; 				//chunksize for the last node [possible exception]
		}

	// Verify structure and reduce
	for (i = 0; i < size; i++) {
		if(myid ==0){
			if (aug_matrix[i][i] == 0) {
			//If the element i,i is a zero, iterate through the rows until finding
			//one whose element is non-zero, then swaps them. If none are found,
			//program exits, as matrix does not have a unique solution.
				for (j = i + 1; j < size ; j++) {
					if (aug_matrix[j][i] != 0) {
					// [j][i] por que no contamos la pivote solo buscamos remplazo q tenga la columna no cero
						swap(aug_matrix[i], aug_matrix[j], size);
						break;
					}
				}
				if (j >= size) {
					cout << endl <<"System does not have a unique solution" << endl;
					exit(1);
				}
			}
		}

	}

	// forwards elimination
	for(i = 0; i < size; i++){
		if(myid == 0){
			//If current pivot coefficient is greater than one, multiply row by multiplicative inverse.
			if (aug_matrix[i][i] != 1) {
				inv = modp.mult_inv[aug_matrix[i][i]];
				for (j = 0; j < size + 1; j++) {     //parallel for
					aug_matrix[i][j] = (aug_matrix[i][j] * inv) % modp.p;
				}
			}

			temp_array_pivot = aug_matrix[i];
			k = 1;
		}

			MPI::COMM_WORLD.Bcast(temp_array_pivot, size + 1, MPI::INT, 0);
			MPI::COMM_WORLD.Bcast(&i, 1, MPI::INT, 0);

		for (j = i + 1; j < size; j++){

			MPI::COMM_WORLD.Bcast(&k, 1, MPI::INT, 0);

			if(myid == 0){
					MPI::COMM_WORLD.Send (aug_matrix[j], size+1, MPI::INT, k, tag);
					MPI::COMM_WORLD.Send (&j, 1, MPI::INT, k, tag); //Intenta guradar la posicion especifica de la fila
			}

			if(myid == k){ // k is going to depend on how many n there are.
				MPI::COMM_WORLD.Recv(temp_array,size + 1, MPI::INT, 0, tag);
				MPI::COMM_WORLD.Recv(&j_row, 1, MPI::INT, 0, tag); //Intenta guradar la posicion especifica de la fila

				eliminate(temp_array_pivot, temp_array, i, size, modp); //el primer parametro es lo que recibas por broadcast a temp_array_pivot

				MPI::COMM_WORLD.Send (&j_row, 1, MPI::INT, 0, tag); //Intenta guradar la posicion especifica de la fila
				MPI::COMM_WORLD.Send (temp_array, size + 1, MPI::INT, 0, tag);


			}
							if(myid == 0){
								k++;
								if(k == numprocs){
									k = 1;
								}
							}


		}

		//MPI::COMM_WORLD.Barrier();  //Test if nonblockin or something

		if(myid==0){
			k = 1;
		}
		for (j = i + 1; j < size; j++){
			//Send de vuelta el temp_array modificado
			//id ==0 reiceve a temp_array
			if(myid == 0){
				MPI::COMM_WORLD.Recv(&j_row, 1, MPI::INT, k, tag); //Intenta guradar la posicion especifica de la fila
				MPI::COMM_WORLD.Recv(aug_matrix[j_row],size + 1, MPI::INT, k, tag); // Cambiar por j_row
			}

			if(myid == 0){
				k++;
				if(k == numprocs){
					k = 1;
				}

			}
		}

		MPI::COMM_WORLD.Barrier();

	}

	MPI::COMM_WORLD.Barrier();

	//backwards elimination
	for (i = size-1; i >= 1; i--) {

		if(myid == 0){
			temp_array_pivot = aug_matrix[i];
			k=1;
		}

		//Broadcast la fila pivot actual aqui
		MPI::COMM_WORLD.Bcast(temp_array_pivot, size + 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&i, 1, MPI::INT, 0);

		for(j = i-1; j >= 0; j--) {

			MPI::COMM_WORLD.Bcast(&k, 1, MPI::INT, 0);

			if(myid == 0){
					MPI::COMM_WORLD.Send (aug_matrix[j], size+1, MPI::INT, k, tag);
					MPI::COMM_WORLD.Send (&j, 1, MPI::INT, k, tag); //Intenta guradar la posicion especifica de la fila
			}

			if(myid == k){ // k is going to depend on how many n there are.
				MPI::COMM_WORLD.Recv(temp_array,size + 1, MPI::INT, 0, tag);
				MPI::COMM_WORLD.Recv(&j_row, 1, MPI::INT, 0, tag); //Intenta guradar la posicion especifica de la fila

				eliminate(temp_array_pivot, temp_array, i, size, modp); //el primer parametro es lo que recibas por broadcast a temp_array_pivot

				MPI::COMM_WORLD.Send (&j_row, 1, MPI::INT, 0, tag); //Intenta guradar la posicion especifica de la fila
				MPI::COMM_WORLD.Send (temp_array, size + 1, MPI::INT, 0, tag);


			}
							if(myid == 0){
								k++;
								if(k == numprocs){
									k = 1;
								}
							}

		}

		if(myid==0){
			k = 1;
		}

		for(j = i-1; j >= 0; j--) {
			//Send de vuelta el temp_array modificado
			//id ==0 reiceve a temp_array
			if(myid == 0){
				MPI::COMM_WORLD.Recv(&j_row, 1, MPI::INT, k, tag); //Intenta guradar la posicion especifica de la fila
				MPI::COMM_WORLD.Recv(aug_matrix[j_row],size + 1, MPI::INT, k, tag); // Cambiar por j_row
			}

			if(myid == 0){
				k++;
				if(k == numprocs){
					k = 1;
				}

			}
		}

		MPI::COMM_WORLD.Barrier();

	}

	MPI::COMM_WORLD.Barrier();

	return;
}


void Matrix::gauss(Space modp) {

	int numprocs, myid, tag = 0;

    numprocs = MPI::COMM_WORLD.Get_size();
    myid = MPI::COMM_WORLD.Get_rank();


	int i, j, k, inv;
	int temp_array[size + 1], *temp_array_pivot = new int[size + 1];
	//bool cont = false;

	for (i = 0; i < size; i++) {
		if(myid ==0){
		if (aug_matrix[i][i] == 0) {
		//If the element i,i is a zero, iterate through the rows until finding
		//one whose element is non-zero, then swaps them. If none are found,
		//program exits, as matrix does not have a unique solution.
			for (j = i + 1; j < size; j++) {
				if (aug_matrix[j][i] != 0) {
					swap(aug_matrix[i], aug_matrix[j], size);
					break;
				}
			}
			if (j >= size) {
				cout << endl <<"System does not have a unique solution" << endl;
				exit(1);
			}
		}
		if (aug_matrix[i][i] != 1) {
		//If pivot coefficient is greater than one, multiply row by multiplica-
		//tive inverse.
			inv = modp.mult_inv[aug_matrix[i][i]];
			for (j = 0; j < size + 1; j++) {     //parallel for
				aug_matrix[i][j] = (aug_matrix[i][j] * inv) % modp.p;
			}
		}

		temp_array_pivot = aug_matrix[i];
	}
		MPI::COMM_WORLD.Bcast(temp_array_pivot, size + 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&i, 1, MPI::INT, 0);

		if(myid == 0){
			k = 1;
		}
		for (j = i + 1; j < size; j++){

			MPI::COMM_WORLD.Bcast(&k, 1, MPI::INT, 0);

		//Reduces rows.
			//id ==0 Send con el row, el size, y el modp
			if(myid == 0){ //
					MPI::COMM_WORLD.Send (aug_matrix[j], size+1, MPI::INT, k, tag);
			}

				//id == j%numprocs Recieve;
			if(myid == k){
				MPI::COMM_WORLD.Recv(temp_array,size + 1, MPI::INT, 0, tag);

				eliminate(temp_array_pivot, temp_array, i, size, modp); //el primer parametro es lo que recibas por broadcast a temp_array_pivot

				MPI::COMM_WORLD.Send (temp_array, size + 1, MPI::INT, 0, tag);


			}

			if(myid == 0){
				k++;
				if(k == numprocs){
					k = 1;
				}
			}

		}
		//MPI::COMM_WORLD.Barrier();

		if(myid == 0){
			k=1;
		}
		for (j = i + 1; j < size; j++){
			//Send de vuelta el temp_array modificado
			//id ==0 reiceve a temp_array
			if(myid == 0){
				MPI::COMM_WORLD.Recv(aug_matrix[j],size + 1, MPI::INT, k, tag);
				for(int o =0; o < size+1; o++){
					cout << aug_matrix[j][o] << " ";
				}
				cout << endl;
			}

			if(myid == 0){
				k++;
				if(k == numprocs){
					k = 1;
				}
				display();
			}
		}

		MPI::COMM_WORLD.Barrier();
	}
	MPI::COMM_WORLD.Barrier();




	for (i = size-1; i >= 1; i--) {

		//Broadcast la fila pivot actual aqui
		temp_array_pivot = aug_matrix[i];
		MPI::COMM_WORLD.Bcast(temp_array_pivot, size + 1, MPI::INT, 0);

		k=1;
		for(j = i-1; j >= 0; j--) {
		/*
		cont = false;
				if(myid == 0){
					if( aug_matrix[j][j] != 0) {
						cont = true;
					}
				}
				*/
		//MPI::COMM_WORLD.Bcast(&cont, 1, MPI::BOOL, 0);
				//      if(cont){

				//id ==0 Send con el row, el indice, el size, y el modp
				if(myid == 0){
					MPI::COMM_WORLD.Send (aug_matrix[j], size + 1, MPI::INT, k, tag);
				}
				//id == j%numprocs Recieve;
				if(myid == k){
					MPI::COMM_WORLD.Recv(temp_array,size + 1, MPI::INT, 0, tag);

					eliminate(temp_array_pivot, temp_array, i, size + 1, modp); //el primer parametro es lo que recibas por broadcast a temp_array_pivot

					MPI::COMM_WORLD.Send (temp_array, size + 1, MPI::INT, 0, tag);
				}
				k++;
				if(k == numprocs){
					k = 1;
				}
			//      }
		}
		//MPI::COMM_WORLD.Barrier();

		k = 1;
		for(j = i-1; j >= 0; j--) {
			/*
			cont = false;
				if(myid == 0){
					if( aug_matrix[j][j] != 0) {
						cont = true;
					}
				}
				*/
			//MPI::COMM_WORLD.Bcast(&cont, 1, MPI::BOOL, 0);
				//       if(cont){

				//Send de vuelta el temp_array modificado
				//id ==0 reiceve a temp_array
				if(myid == 0){
				MPI::COMM_WORLD.Recv(aug_matrix[j],size + 1, MPI::INT, k, tag);
				for(int o =0; o < size+1; o++){
						cout << aug_matrix[j][o] << " ";
					}
					cout << endl;
				}
				k++;
				if(k == numprocs){
					k = 1;
				}
			//       }
			if(myid == 0){
				display();
			cout << endl << endl << endl;
			}
		}
		MPI::COMM_WORLD.Barrier();
	}

	return;
}

