/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#include "teaser/graph.h"
#include "pmc/pmc.h"
#include <metis.h>
#include <mpi.h>
#include <vector>
//#include <boost/mpi/datatype.hpp>
using namespace std;

vector<int> teaser::MaxCliqueSolver::findMaxClique(teaser::Graph graph) {

  // Handle deprecated field
  if (!params_.solve_exactly) {
    params_.solver_mode = CLIQUE_SOLVER_MODE::PMC_HEU;
  }

  // Create a PMC graph from the TEASER graph
  vector<int> edges;
  vector<long long> vertices;// row ptr

  vertices.push_back(edges.size()); // TODO:confirm edges.size() at this point is zero

  // can't be parallellised because  vertices are ordered by weight : refer to build/pmc-src/README.md
  const auto all_vertices = graph.getVertices();
  for (const auto& i : all_vertices) {
    const auto& c_edges = graph.getEdges(i);
    edges.insert(edges.end(), c_edges.begin(), c_edges.end());
    vertices.push_back(edges.size()); 
  }

  // TODO::GATHER local vertices

  // Use PMC to calculate
  pmc::pmc_graph G(vertices, edges);
   int numproc,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  cout << "rank search dense****** : " << rank << endl;
  vector<int> local_rowptr;
  vector<int> local_colind;
  vector<int> partition_nedges(numproc);
  vector<int>  partition_nvertices(numproc);
 
  //if (rank == 0){

        idx_t num_of_vertices =  G.vertices.size();
        cout << "vertex size****** : " << num_of_vertices<< endl;
        idx_t num_of_edges =  G.edges.size();
        cout << "edges size****** : " << num_of_edges<< endl;
       // graph to a CSR (Compressed Sparse Row) format
        idx_t *rowptr = new idx_t [num_of_vertices ] ;
        idx_t *colind = new idx_t [num_of_edges] ;

        cout << "Saving rowpointers and column indices to file" << endl;
        ofstream myrpfile;
        myrpfile.open ("rowptr.txt");
        if( !myrpfile ) { // file couldn't be opened
        cerr << "Error: rp file could not be opened" << endl;
        exit(1);
        }
        ofstream mycolindfile;
        mycolindfile.open("colind.txt");
        if( !mycolindfile ) { // file couldn't be opened
        cerr << "Error: colind file could not be opened" << endl;
        exit(1);
        }
        ofstream mypartrpfile;
        ofstream mypartcolindfile;

        for (int j = 0; j< num_of_edges ; j++) {
           colind[j] = (idx_t) G.edges[j];
        }

        for (int m = 0; m< num_of_vertices ; m++){
            rowptr[m] = (idx_t) G.vertices[m] ;
        }

        for (int i = 0; i< num_of_vertices ; i++) {     
        myrpfile << rowptr[i] << endl ;
        for (int c = rowptr[i] ; c < rowptr[i+1]; c++){
            mycolindfile << colind[c] << " " ;     
        }
        mycolindfile << endl;
        }
        myrpfile.close();
        mycolindfile.close();
        cout << "files closed " << endl;
        cout << "before partition " << endl;

        // Set up the options for the partitioning
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options); 
        options[METIS_OPTION_DBGLVL] = 255; // adjust debugging level of metis
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
        options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
        options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
        options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
        options[METIS_OPTION_NCUTS] = 1;
        options[METIS_OPTION_NITER] = 10;
        //TODO : check options[METIS OPTION NUMBERING]
        idx_t nparts = 8; 
        cout << "nparts : " << nparts  << endl;
        cout << "n vertices : " << num_of_vertices  << endl;
        cout << "n edges : " << num_of_edges  << endl;
        idx_t objval ;
        idx_t constraints = 1;
        idx_t *part = new idx_t [num_of_vertices];
        num_of_vertices--;
        cout << "partitioning...  " << endl;
        int status =  METIS_PartGraphKway(&num_of_vertices,&constraints  /*balance constraints*/, rowptr, colind,NULL, NULL, NULL, &nparts,
        NULL,NULL,options, &objval/* output NO OF EDGES OF PART*/, part /*shows which vertex has been assigned to which part*/);
        cout << "objval : " << objval << endl;
        cout << "partition done" << endl;

        /*
        TODO:
          1. create 2D array for Partition column indices [number of parts] [adjacecyvals for that partition]
          2. create 2D array for Partition rowptrs [number of parts] [ptrs for that partition ]
          3. distribute the 1D arrays to the respective processes
        
        */
      mypartcolindfile.open("partcolind.txt");
      if( !mypartcolindfile) { // file couldn't be opened
        cerr << "Error: partition colind file could not be opened" << endl;
        exit(1);
      }
      mypartrpfile.open ("partrowptr.txt");
      if( !mypartrpfile ) { // file couldn't be opened
        cerr << "Error:partition  rp file could not be opened" << endl;
        exit(1);
      }
    
        
    // store partitions in 2D arrays of row ptr and column ids
   
      vector<vector<idx_t>> partition_rowptr(nparts);
      vector<vector<idx_t>> partition_colidx(nparts); 

      vector<int> maps(num_of_vertices);
      vector<int> cnt(nparts,0);
      vector<int> nz(nparts,0);
      vector<int> rpcounter(nparts,1);
  
      int sz = num_of_vertices;
      //initialise to 0
     for (int p = 0; p < nparts; p++){
        partition_nedges[p] = 0;
        partition_nvertices[p] = 0;
      }



     for (int i = 0; i < sz ; i++) {
        partition_nvertices[part[i]]++;
        maps[i]= cnt[part[i]]++;
        for (int j = rowptr[i]; j < rowptr[i + 1]; j++) {
             if (part[colind[j]]== part[i]){
                partition_nedges[part[i]]++;
                }
        }

      }

    for (int p = 0; p < nparts; p++) {
        partition_rowptr[p].resize(partition_nvertices[p]+1);
        partition_colidx[p].resize(partition_nedges[p]);
        partition_rowptr[p][0] = 0;
    }
    cout << "OK 1" << endl;
    
    for (int i = 0; i < sz ; i++) {
        int p = part[i]; // Get partition of current vertex
        for (int j = rowptr[i]; j < rowptr[i + 1]; j++) {
            int neighbor = colind[j];

            if (part[neighbor]== p){
                partition_colidx[p][nz[p]] = maps[neighbor];
                nz[p]++;
                }
   
        }

        partition_rowptr[p][rpcounter[p]++] = nz[p];

    }
     cout << "OK 2" << endl;

     cout << "Partition csr structure created" << endl;
    // save to file
    for (int r = 0 ; r < nparts ; r++){
      mypartrpfile << "part "<< r << endl;  
      for (int w = 0 ; w < partition_rowptr[r].size() ; w++){
          mypartrpfile << partition_rowptr[r][w] << " ";
        }
       mypartrpfile << endl;
       mypartrpfile << endl;
     }

    for (int r = 0 ; r < nparts ; r++){
      mypartcolindfile << "part "<< r << endl;  
      for (int w = 0 ; w < partition_colidx[r].size() ; w++){
          mypartcolindfile<< partition_colidx[r][w] << " ";
        }
       mypartcolindfile << endl;
       mypartcolindfile << endl;
     }

      mypartcolindfile.close();
      mypartrpfile.close();
    
    
    // vector to represent max clique
    vector<vector<int>> C(nparts);

    for (int r = 0 ; r< nparts;r++){
     vector<long long> rplong(begin(partition_rowptr[r]), end(partition_rowptr[r]));
     pmc::pmc_graph local_G(rplong,partition_colidx[r]);
     pmc::input local_in;
      local_in.algorithm = 0;
      local_in.threads = 8;
      local_in.experiment = 0;
      local_in.lb = 0;
      local_in.ub = 0;
      local_in.param_ub = 0;
      local_in.adj_limit = 20000; 
      local_in.time_limit = params_.time_limit;
      local_in.remove_time = 4;
      local_in.graph_stats = false;
      local_in.verbose = false;
      local_in.help = false;
      local_in.MCE = false;
      local_in.decreasing_order = false;
      local_in.heu_strat = "kcore";
      local_in.vertex_search_order = "deg";


    // upper-bound of max clique
    local_G.compute_cores();
    auto loc_max_core = local_G.get_max_core();
    //auto mdt = boost::mpi::get_mpi_datatype(max_core);
    TEASER_DEBUG_INFO_MSG("Max core number: " << loc_max_core);
    TEASER_DEBUG_INFO_MSG("Num vertices: " << local_G.vertices.size() - 1); // ?? recheck
    //}
    //MPI_Bcast(&max_core,1,mdt,0,MPI_COMM_WORLD);


  // check for k-core heuristic threshold
  // check whether threshold equals 1 to short circuit the comparison
  if (params_.solver_mode == CLIQUE_SOLVER_MODE::KCORE_HEU &&
      params_.kcore_heuristic_threshold != 1 &&
      loc_max_core > static_cast<int>(params_.kcore_heuristic_threshold *
                                  static_cast<double>(local_G.vertices.size()))) {
    TEASER_DEBUG_INFO_MSG("Using K-core heuristic finder.");
    // remove all nodes with core number less than max core number
    // k_cores is a vector saving the core number of each vertex
    auto k_cores = local_G.get_kcores();
    for (int i = 1; i < k_cores->size(); ++i) {
      // Note: k_core has size equals to num vertices + 1
      if ((*k_cores)[i] >= loc_max_core) {
        C[r].push_back(i-1);
      }
    }
    return C[r];
  }

  if (local_in.ub == 0) {
    local_in.ub = loc_max_core + 1;
  }

  // lower-bound of max clique
  if (local_in.lb == 0 && local_in.heu_strat != "0") { // skip if given as input
    pmc::pmc_heu maxclique(local_G, local_in);
    local_in.lb = maxclique.search(local_G, C[r]);
  }

  assert(local_in.lb != 0);
  if (local_in.lb == 0) {
    // This means that max clique has a size of one
    TEASER_DEBUG_ERROR_MSG("Max clique lower bound equals to zero. Abort.");
    return C[r];
  }

  if (local_in.lb == local_in.ub) {
    return C[r];
  }
  if (params_.solver_mode == CLIQUE_SOLVER_MODE::PMC_EXACT) {

      int mcval;
      if (local_G.num_vertices() < local_in.adj_limit) {
          local_G.create_adj();
          //cout << "OK b" << endl;
          pmc::pmcx_maxclique finder(local_G, local_in);
          //cout << "OKkkk" << endl;
          mcval = finder.search_dense(local_G, C[r]);
          //cout << "OK c" << endl;
        } else {
          pmc::pmcx_maxclique finder(local_G, local_in);
          mcval = finder.search(local_G, C[r]);
        }
        cout << "mcval " << r << " " << mcval << endl;
    }

  }

  return C[0];
  // // Prepare PMC input
  // // TODO: Incorporate this to the constructor
  // pmc::input in;
  // in.algorithm = 0;
  // in.threads = 8;
  // in.experiment = 0;
  // in.lb = 0;
  // in.ub = 0;
  // in.param_ub = 0;
  // in.adj_limit = 20000; 
  // in.time_limit = params_.time_limit;
  // in.remove_time = 4;
  // in.graph_stats = false;
  // in.verbose = false;
  // in.help = false;
  // in.MCE = false;
  // in.decreasing_order = false;
  // in.heu_strat = "kcore";
  // in.vertex_search_order = "deg";

  // // vector to represent max clique
  // vector<int> C;

  //   // upper-bound of max clique
  //   G.compute_cores();
  //   auto max_core = G.get_max_core();
  //   //auto mdt = boost::mpi::get_mpi_datatype(max_core);
  //   TEASER_DEBUG_INFO_MSG("Max core number: " << max_core);
  //   TEASER_DEBUG_INFO_MSG("Num vertices: " << vertices.size());
  //   //}
  //   //MPI_Bcast(&max_core,1,mdt,0,MPI_COMM_WORLD);


  // // check for k-core heuristic threshold
  // // check whether threshold equals 1 to short circuit the comparison
  // if (params_.solver_mode == CLIQUE_SOLVER_MODE::KCORE_HEU &&
  //     params_.kcore_heuristic_threshold != 1 &&
  //     max_core > static_cast<int>(params_.kcore_heuristic_threshold *
  //                                 static_cast<double>(all_vertices.size()))) {
  //   TEASER_DEBUG_INFO_MSG("Using K-core heuristic finder.");
  //   // remove all nodes with core number less than max core number
  //   // k_cores is a vector saving the core number of each vertex
  //   auto k_cores = G.get_kcores();
  //   for (int i = 1; i < k_cores->size(); ++i) {
  //     // Note: k_core has size equals to num vertices + 1
  //     if ((*k_cores)[i] >= max_core) {
  //       C.push_back(i-1);
  //     }
  //   }
  //   return C;
  // }

  // if (in.ub == 0) {
  //   in.ub = max_core + 1;
  // }

  // // lower-bound of max clique
  // if (in.lb == 0 && in.heu_strat != "0") { // skip if given as input
  //   pmc::pmc_heu maxclique(G, in);
  //   in.lb = maxclique.search(G, C);
  // }

  // assert(in.lb != 0);
  // if (in.lb == 0) {
  //   // This means that max clique has a size of one
  //   TEASER_DEBUG_ERROR_MSG("Max clique lower bound equals to zero. Abort.");
  //   return C;
  // }

  // if (in.lb == in.ub) {
  //   return C;
  // }

  // // Optional exact max clique finding
  // if (params_.solver_mode == CLIQUE_SOLVER_MODE::PMC_EXACT) {
  //   // The following methods are used:
  //   // 1. k-core pruning
  //   // 2. neigh-core pruning/ordering
  //   // 3. dynamic coloring bounds/sort
  //   // see the original PMC paper and implementation for details:
  //   // R. A. Rossi, D. F. Gleich, and A. H. Gebremedhin, “Parallel Maximum Clique Algorithms with
  //   // Applications to Network Analysis,” SIAM J. Sci. Comput., vol. 37, no. 5, pp. C589–C616, Jan.
  //   // 2015.

  //   //********************************************

 
  //   //***********************************************
  //   for (int r = 0 ; r< nparts;r++){
  //     vector<long long> rplong(begin(partition_rowptr[r]), end(partition_rowptr[r]));
  //     pmc::pmc_graph local_G(rplong,partition_colidx[r]);
  //     local_G.compute_cores();
  //     // vector to represent max clique
  //     int mcval;
  //     if (local_G.num_vertices() < in.adj_limit) {
  //       local_G.create_adj();
  //       //cout << "OK b" << endl;
  //       pmc::pmcx_maxclique finder(local_G, in);
  //       //cout << "OKkkk" << endl;
  //       mcval = finder.search_dense(local_G, C);
  //       //cout << "OK c" << endl;
  //     } else {
  //       pmc::pmcx_maxclique finder(local_G, in);
  //       mcval = finder.search(local_G, C);
  //     }
  //     cout << "mcval " << r << " " << mcval << endl;
  //   }
  
  // }

  

}