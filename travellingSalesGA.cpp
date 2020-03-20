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
#define CROSSBREED 16
#define XBREEDSHR 16
#define MATINGPOOL 256
#define MUTATIONRT 0.02
#define POPULATION 1024

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
      // cout << "while(getline(fs, line))\n";
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
    // cout << "\n\n";
    cout << "----------Error in swap()!----------\n";
    // for(int i = 0; i < seq.size(); ++i){
    //   cout << seq[i].first << ", " << seq[i].second << "\n";
    // }
    // cout << "\ncoord = (" << coord.first << ", " << coord.second << ")\n"; //TODO: cut this out
    // cout << "vector index = " << idx <<"\n";
    // cout << "seq[idx] = (" << seq[idx].first << ", " << seq[idx].second << ")\n";
    // cout << "seq[it] = (" << it->first << ", " << it->second << ")\n";
  }
} /* swap() */

bool selection(vector<vector<pair<int,int>>> &seqs, int srvvrs){
  bool converged = false;
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

  /* while the candidate list has not reached the desired size */
  while(sz > srvvrs){
    /* implement a stochastic acceptance selection */
    int idx = rand()%sz;
    double rand0to1 = ((double)rand() / (double)RAND_MAX);

    /* determine whether to discard this sequence */
    if(((dist[idx]-min_dist)/(max_dist-min_dist)) > rand0to1 ||
       max_dist == min_dist){
      if(max_dist == min_dist){
        converged = true;
      }
      dist_sum -= dist[idx]; /* decrement sum of distances */
      sz--; /* decrement size */
      dist.erase(dist.begin()+idx); /* delete element */
      seqs.erase(seqs.begin()+idx); /* delete element */
    } else {
      /* check for complete convergence by recalculating max_dist */
      max_dist = 0.0;
      for(int i = 0; i < sz; ++i){
        if(dist[i] > max_dist){
          max_dist = dist[i];
        }
      }
    } /* end else not discarding this sequence */
  } /* end while size is greater than the number of survivors */

  return converged;
} /* selection() */

vector<pair<int,int>> pmx(vector<pair<int,int>> &seq1,
                          vector<pair<int,int>> seq2){
  /* Partially Matched Crossover (PMX) function returns a child
     sequence from the two parent sequences.
     It assumes that the sequences are of the same size            */
// cout << seq1.size() << " = " << seq2.size() << "\n";
  vector<pair<int,int>>::size_type start, end;/* indexes for gene */
  start = rand()%(seq2.size()-1);/* pick a random index to start */
  end = rand()%(seq2.size()-(start+1))+start+1;/*and a ending index*/

  /* perform crossover */
  for(int i = start; i<=end; ++i){
    // cout << "pmx()\n";
    swap(seq2, seq1[i], i);
  }
  return seq2;
} /* pmx() */

void mutate(vector<pair<int,int>> &seq){
  /* This mutate function will mutate a sequence by 2 methods:
     1. shifting the sequence forward (with wraparound)
     2. swapping values in the sequence                            */

  /* shift sequence (with wraparaound) */
  pair<int,int> tmp;
  for(int i = 0; i < (rand()%5+1); ++i){
    tmp = seq.front();
    seq.erase(seq.begin());
    seq.push_back(tmp);
  }

  /* swap values in sequence */
  for(int i = 0; i < (rand()%5+1); ++i){
    int idx1 = rand()%seq.size();
    int idx2 = rand()%seq.size();
    // cout << "mutate()\n";
    swap(seq, seq[idx1], idx2);
  }
} /* mutate() */

void reproduce(vector<vector<pair<int,int>>> &seqs, int numchild){
  /* This reproduce function will add children sequences to the
     current population of sequences.                              */
  int newchildren = 0;
  while(newchildren < numchild){
    // cout<<"while(newchildren < numchild)\n";
    int prnt1 = rand()%seqs.size(); /* randomly pick a parent */
    int prnt2 = rand()%seqs.size(); /* prnt1==prnt2 is acceptable */
    seqs.push_back(pmx(seqs[prnt1], seqs[prnt2]));
    newchildren++;

    /* perform mutation TODO*/
    if(MUTATIONRT > ((double)rand() / (double)RAND_MAX)){
      mutate(seqs.back());
    }
  }
} /* reproduce() */

int xbreed(vector<vector<pair<int,int>>> &seqs, int rank, int np){
  /* This crossbreed function serves to mix up genes between the
     several processors that are all running on this MPI
     communicator.                                                 */
  int sz = seqs.size(); /* total number of sequences */
  int xbred = 0; /* number of crossbreeding events */
  int shareseqs[XBREEDSHR*2]; /* int array for sharing sequences */

  /* determine sequences to share with others */
  vector<vector<pair<int,int>>> subseqs = seqs;
// cout << subseqs.size() << "\n";
  selection(subseqs, XBREEDSHR);
// cout << subseqs.size() << "\n";
  for(int i = 0; i < XBREEDSHR; ++i){/* for each sequence to share */
    vector<pair<int,int>> seq; /* temporary sequence */

    /* convert sequence into an array of integers {X1, Y1, X2,...} */
    for(int j = 0; j < subseqs[j].size(); ++j){/* for each pair */
      shareseqs[j*2] = subseqs[i][j].first;/* x-coordinate */
      shareseqs[j*2+1] = subseqs[i][j].second;/* y-coordinate */
    }

    /* share sequence with next higher rank processor (wraparound) */
    MPI_Send(shareseqs, XBREEDSHR*2, MPI_INT, (rank+1)%np, 0, MCW);
    /* receive sequence from next lower rank processor (wraparound)*/
    MPI_Recv(shareseqs, XBREEDSHR*2, MPI_INT, (rank+np-1)%np, 0, MCW,
        MPI_STATUS_IGNORE);

    /* convert received ints into a sequence */
    for(int j = 0; j < subseqs[j].size(); ++j){
      seq.push_back(make_pair(shareseqs[j*2], shareseqs[j*2+1]));
      // cout << shareseqs[j*2] << ", " << shareseqs[j*2+1] << "\n";
    }
    seqs.push_back(seq); /* add this sequence to the big list */
    ++xbred; /* increment cross breeding counter */
  } /* end for each sequence to share */

  /* get seqs back to normal size */
  selection(seqs, sz);

  return xbred;
} /* xbreed() */

bool generation(vector<vector<pair<int,int>>> &seqs, int numparnts){
  /* This generation function will create a new generation with
     surviviors of the previous generation and their surviving 
     children sequences.
     numparnts must be even.                                       */
  int sz = seqs.size();
  if(numparnts%2){
    cout << "Error in generation(): uneven number of parents\n";
    return true;
  } else if(sz < POPULATION || numparnts > POPULATION){
    cout << "Error in generation(): invalid size of populace.\n";
    return true;
  }

  /* select parents (reduce the number of sequences) */
  selection(seqs, numparnts);

  /* create children from parents (increase number beyond sz) */
  reproduce(seqs, sz);

  /* determine survivors (cuts number back down to sz) */
  return selection(seqs, sz);
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

  /* create candidates */
  for(int i = 0; i < POPULATION; ++i){
    sequences.push_back(sequence);
  }

  /* shuffle sequences */
  for(int i = 0; i < POPULATION; ++i){
    random_shuffle(sequences[i].begin(), sequences[i].end());
  }

  double max_fit = 100000000000000.0;
  int running = 128;
  int iters = 1;
  int xbred = 0;

  while(running){
  
    /* fitness measure */
    for(int i=0; i < sequences.size(); ++i){
      double fit = total_dist(sequences[i]);
      if(fit < max_fit){max_fit = fit;}
    }

    /* report current status */
    cout << "distance(" << iters++ << ") = " << max_fit << ";\n";

    /* iterate a generation */
    bool converged = generation(sequences, MATINGPOOL);
    // if(converged){
    //   cout << rank << "'s sequence has converged!\n";
    // }
    /* swap some data around (crossbreed) */
    if(!(iters%CROSSBREED) || converged){ /* if time to crossbreed */
cout << sequences.size() << " = ";
      xbred += xbreed(sequences, rank, np);
cout << sequences.size() << "\n";
    }

    /* check for stopping condition */
    running--;
  } /* end while(running) */

  /* terminate program */
  MPI_Finalize();
  return 0;
} /* main() */
