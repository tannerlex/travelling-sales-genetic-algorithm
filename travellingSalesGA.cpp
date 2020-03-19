#include <algorithm>
#include <fstream>
#include <iostream>
// #include <map>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <time.h> /* used to seed rand with current time */
#include <vector>

using namespace std;

// #define INFILE "locations.txt"
#define INFILE "locations.txt"
#define MCW MPI_COMM_WORLD
#define MATINGPOOL 256
#define POPULATION 1024
#define CROSSBREED 16

bool read_locations(string &filename, vector<pair<int,int>> &data){
/* file must be in the form:
XXXXXXXX YYYYY
 XXXXXXX YYYYYY
etc. where each line consists of 2 whitespace separated integers
(a comma separation is also acceptable as stoi() will ignore it) */

  string line; /* used to store a single line in the file */
  bool success; /* return flag indicating success of file read */
  ifstream fs(filename);/* stream for reading file */
  if(fs.is_open()){ /* input file was opened */
    while(getline(fs, line)){
      stringstream ss(line); /* stream for each line */
      string x,y;

      /* read and record each pair of ints */
      if(ss >> x >> y){
        /* add each coordinate pair to vector */
        data.push_back(make_pair(stoi(x),stoi(y)));
      }
    } /* end while streaming file */
    success = true;
    fs.close();
  } else { /* input file was not opened */
    cout << "File: \"" << filename << "\" not found.\n";
    success = false;
  }
  return success;
} /* read_locations() */

void print_seq(vector<pair<int,int>> &seq){
  for(auto & coord : seq){
    cout << "(" << coord.first << ", " << coord.second << ")\n";
  }
} /* print_seq() */

double total_dist(vector<pair<int,int>> &seq){
  double length = 0.0; /* running total of sequence distance */

  /* for each pair of coordinates */
  for(vector<pair<int,int>>::size_type i=0; i!=(seq.size()-1); ++i){
    /* add the distance between each pair of coordinates to total */
    length += sqrt( /* distance is Pythagorean theorem */
        pow((double)(seq[i].first - seq[i+1].first), 2.0) +
        pow((double)(seq[i].second - seq[i+1].second), 2.0));
  }
  return length;
} /* total_dist() */

void swap(vector<pair<int,int>> &seq, pair<int,int> &coord, int idx){
  /* Swap the coordinate at the specified index with the specified
     coordinate wherever it may be within the sequence             */
  auto it = find(seq.begin(), seq.end(), coord);
  if (it != seq.end()){
    std::swap(seq[idx], *it);
  } else {
    cout << "----------Error in swap()!----------\n";
    // cout << "coord = (" << coord.first << ", " << coord.second << ")\n"; //TODO: cut this out
    // cout << "vector index = " << idx <<"\n";
    // cout << "seq[idx] = (" << seq[idx].first << ", " << seq[idx].second << ")\n";
    // cout << "seq[it] = (" << it->first << ", " << it->second << ")\n";
  }
} /* swap() */

void selection(vector<vector<pair<int,int>>> &seqs, int srvvrs){
  int sz = seqs.size();
  double dist_sum = 0.0;
  vector<double> dist;
  double max_dist = 0.0;
  double min_dist = 100000000000000.0;
  for(int i=0; i < sz; ++i){
    dist.push_back(total_dist(seqs[i]));
    dist_sum += dist[i];
    if(dist[i]>max_dist){
      max_dist = dist[i];
    }
    if(dist[i]<min_dist){
      min_dist = dist[i];
    }
  }
  // double av = dist_sum/sz;
  // cout << "BEFORE\n";
  // cout << "total:   " << dist_sum << "\n";
  // cout << "size:    " << sz << "\n";
  // cout << "average: " << av << "\n\n";
  // int morethan = 0;
  // int lessthan = 0;
  while(sz > srvvrs){
    int idx = rand()%sz;
    /* implement a stochastic acceptance selection */
    double rand0to1 = ((double)rand() / (double)RAND_MAX);

    /* determine whether to discard this sequence */
    if(((dist[idx]-min_dist)/(max_dist-min_dist))>rand0to1){
    // if((dist[idx]/max_dist)>rand0to1){
      // if(dist[idx]>av){ /* testing the selection algorithm TODO: remove*/
      //   morethan++;
      // } else {lessthan++;}
      // cout << ((dist[idx]-min_dist)/(max_dist-min_dist)) << " " << dist[idx] << "\n";
      dist_sum -= dist[idx]; /* decrement sum of distances */
      sz--; /* decrement size */
      dist.erase(dist.begin()+idx);
      seqs.erase(seqs.begin()+idx);
    }
  }
  // cout << morethan << " sequences removed had higher than average distance.\n";
  // cout << lessthan << " sequences removed had less than average distance.\n";
  // cout<<endl;
  // for(int i =0; i<sz; ++i){
  //   cout << dist[i] << " ";
  // }
  // cout << "\nAFTER\n";
  // cout << "total:   " << dist_sum << "\n";
  // cout << "size:    " << sz << "\n";
  // cout << "average: " << dist_sum/sz << "\n";
} /* selection() */

vector<pair<int,int>> pmx(vector<pair<int,int>> &seq1,
                          vector<pair<int,int>> seq2){
  /* Partially Matched Crossover (PMX) function returns a child
     sequence from the two parent sequences.
     It assumes that the sequences are of the same size            */

  vector<pair<int,int>>::size_type start, end;/* indexes for gene */
  start = rand()%(seq2.size()-1);/* pick a random index to start */
  end = rand()%(seq2.size()-(start+1))+start+1;/*and a ending index*/

  /* perform crossover */
  for(int i = start; i<=end; ++i){
    swap(seq2, seq1[i], i);
  }
  return seq2;
} /* pmx() */

void mutate(vector<pair<int,int>> &seq){
  // cout<<seq.size()<<endl;
  pair<int,int> tmp;
  for(int i = 0; i < (rand()%5+1); ++i){
    tmp = seq.front();
    seq.erase(seq.begin());
    seq.push_back(tmp);
  }
  // cout<<seq.size()<<endl;
} /* mutate() */

void reproduce(vector<vector<pair<int,int>>> &seqs, int numchild){
  /* This reproduce function will add children sequences to the
     current population of sequences.                              */
  int newchildren = 0;
  while(newchildren < numchild){
    int prnt1 = rand()%seqs.size(); /* randomly pick a parent */
    int prnt2 = rand()%seqs.size(); /* prnt1==prnt2 is acceptable */
    seqs.push_back(pmx(seqs[prnt1], seqs[prnt2]));
    newchildren++;

    /* perform mutation TODO*/
    if(.05 > ((double)rand() / (double)RAND_MAX)){
      mutate(seqs.back());
    }
  }
} /* reproduce() */

void generation(vector<vector<pair<int,int>>> &seqs, int numparnts){
  /* This generation function will create a new generation with
     surviviors of the previous generation and their surviving 
     children sequences.
     numparnts must be even.                                       */
  int sz = seqs.size();
  if(numparnts%2){
    cout << "Error in generation(): uneven number of parents\n";
    return;
  } else if(sz < POPULATION || numparnts > POPULATION){
    cout << "Error in generation(): invalid size of populace.\n";
    return;
  }

  /* select parents (reduce the number of sequences) */
  selection(seqs, numparnts);

  /* create children from parents (increase number beyond sz) */
  reproduce(seqs, sz);

  /* determine survivors (cuts number back down to sz) */
  selection(seqs, sz);
} /* generation() */

int main(int argc, char **argv){
  /* initialize MPI variables */
  int rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MCW, &rank); 
  MPI_Comm_size(MCW, &np);

  /* initialize algorithm variables */
  string infile = INFILE;
  vector<vector<pair<int, int>>> sequences;
  vector<pair<int,int>> sequence;
  srand(time(0)+rank); /* seed rand with time */

  /* read input file */
  if(!read_locations(infile, sequence)){ 
    /* if file not read terminate program execution */
    MPI_Finalize();
    return 0;
  }

  // /* test that input file is read properly TODO: remove */
  // for(int i=0; i<np; ++i){
  //   MPI_Barrier(MCW);
  //   for(int j=0; j<POPULATION; ++j){
  //     if(i==rank){
  //       print_seq(sequences[j]);
  //       cout << total_dist(sequences[j]) << endl;
  //     }
  //   }
  // }

  /* create candidates */
  for(int i = 0; i < POPULATION; ++i){
    sequences.push_back(sequence);
  }

  /* shuffle sequences */
  for(int i = 0; i < POPULATION; ++i){
    random_shuffle(sequences[i].begin(), sequences[i].end());
  }

  // /* test the PMX function (after shuffle would be good) todo: remove*/
  // sequences.push_back(pmx(sequences[0],sequences[1]));
  // for(int i=0; i<np; ++i){
  //   MPI_Barrier(MCW);
  //   if(i==rank){ 
  //     cout << rank << endl;
  //     for(int j=0; j < sequences.size(); ++j){
  //       print_seq(sequences[j]);
  //       cout << total_dist(sequences[j]) << endl;
  //     }
  //   }
  // }

  /* fitness contest */
  double max_fit = 100000000000000.0;
  vector<double> fitness;
  for(int i = 0; i < sequences.size(); ++i){
    double fit = total_dist(sequences[i]); /* measure fitness */
    fitness.push_back(fit); /* record the fitness of this element */
    if(fit < max_fit){max_fit = fit;} /* keep track of the fittest */

  } /* end for each sequence */

  int running = 32;
  int iters = 1;
  cout << "distance(" << iters++ << ") = " << max_fit << ";\n";
  while(running){
    /* iterate a generation */
    generation(sequences, MATINGPOOL);

    /* swap some data around (crossbreed) */
    if(!(iters%CROSSBREED)){ /* if it's time to crossbreed */
      // xbreed(sequences);
    }

    /* fitness competition */
    for(int i=0; i < sequences.size(); ++i){
      double fit = total_dist(sequences[i]);
      if(fit < max_fit){max_fit = fit;}
    }
    cout << "distance(" << iters++ << ") = " << max_fit << ";\n";

    /* check for stopping condition */
    running--;
  } /* end while(running) */

  /* terminate program */
  MPI_Finalize();
  return 0;
} /* main() */
