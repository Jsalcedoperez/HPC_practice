#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <vector>

namespace defi{

    int fsc_size, new_size;
    std::vector<int> local_sites,sites_per_rank(4);
    int ierr, myrank, tot_ranks;
    int a, b, n, neighbor;
    std::vector<int> target(5), min_work(4);
    int total_sites, remainder, min_target;
    std::vector<MPI_Request> requests;
    
    enum Neighbor {
        LEFT  = 0,
        RIGHT = 1,
        NUM_NEIGHBORS
    };

}


void setup(const std::vector<int>&);

void create_sites();

void rebalance(const std::vector<int>&);

int main(int argc, char *argv[])

{

    using namespace defi;
    int ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &defi::myrank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &defi::tot_ranks);

    std::vector<int> fission_sites = {0, 5, 10, 18, 25};

    setup(fission_sites);
 
    create_sites();
    
    int tot_converged, converge = 0;
    int m = 3;
    if (myrank == m)

    {

    std::cout << "(before) from node "<< myrank <<", this is my (a,b) pair (" << a << ", " << b << ")" << std::endl;
    
    for (int i = 0; i < local_sites.size() ; ++i)
    {

    std::cout << "(before), from node "<<myrank << ", ls[" << i << "] = " << local_sites[i] << std::endl;

    }
  

    }
    while (tot_converged < tot_ranks)
    {

        rebalance(fission_sites);
        fsc_size = b - a + 1;

        if (myrank == m)

        {

        std::cout << "(after) from node "<< myrank <<", this is my (a,b) pair (" << a << ", " << b << ")" << std::endl;
        
        for (int i = 0; i < local_sites.size() ; ++i)
        {
    
        std::cout << "(after), from node "<<myrank << ", ls[" << i << "] = " << local_sites[i] \
        << " and fsc " << fsc_size << std::endl;
        //std::cout << "(after), converge "<< tot_converged << std::endl;

        }
      

        }
        if (fsc_size == min_work[myrank] )
        {

            converge = 1;

        }

         MPI_Allreduce(&converge, &tot_converged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    }


    ierr = MPI_Finalize();

    return 0;
}

void setup(const std::vector<int>& fission_sites)
{
    using namespace defi;
    fsc_size = fission_sites.size();
    std::fill(sites_per_rank.begin(),sites_per_rank.end(),0);
    std::fill(min_work.begin(),min_work.end(),0);
    std::fill(target.begin(),target.end(),0);
    a = fission_sites[myrank];
    b = fission_sites[myrank + 1] -1;

    total_sites = fission_sites[fsc_size-1];
    remainder = total_sites % tot_ranks;
    min_target = total_sites / tot_ranks;
    int ibank = 0;
    for (int irank = 0 ; irank < tot_ranks ; ++irank)
    {
        int imin_work = irank < remainder ? min_target + 1 : min_target;
        ibank += imin_work;
        min_work[irank] = imin_work;
        target[irank + 1] = ibank;

    }

}

void create_sites()
{

    //set up sites
    using namespace defi;
    if (local_sites.empty())
    {

        for (int i = a; i <= b; ++i)
        {
            local_sites.push_back(i);

        }
    }

    
    sites_per_rank[myrank] = b-a+1;
    MPI_Allreduce(MPI_IN_PLACE, sites_per_rank.data(), sites_per_rank.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
    if (sites_per_rank[myrank] != local_sites.size())
    {
       
        std::cout << "sites' size and fsc_size don't match" << std::endl;

    }  

   
}
void rebalance(const std::vector<int>& fission_sites)

{
    int org_size;
    std::cout << "rebalance..." << std::endl;
    using namespace defi;
    //send left
    if (a < target[myrank])
    {

        if (myrank != 0){
        n = std::min((target[myrank]-a), sites_per_rank[myrank]);
        neighbor = myrank - 1;
        requests.emplace_back();
        MPI_Isend(&local_sites[0], n, MPI_INT, neighbor, myrank, MPI_COMM_WORLD, &requests.back());
        a += n;
        sites_per_rank[myrank] -= n;
        }
    }

    //receive from right
    else if (b < (target[myrank + 1]-1))
    {
        if (myrank != (tot_ranks-1)){
  
      
        n = std::min((target[myrank+1]-1-b),sites_per_rank[myrank+1]);
        neighbor = myrank + 1;
        org_size = local_sites.size();
        local_sites.resize(org_size+n);
        requests.emplace_back();
        MPI_Irecv(&local_sites[org_size], n, MPI_INT, neighbor, neighbor, MPI_COMM_WORLD, &requests.back());
        b +=  n;
        sites_per_rank[myrank] += n;

        }

    }

    //receive from left
    if (a > target[myrank])
    {
        n = std::min((a-target[myrank]), sites_per_rank[myrank-1]);
        org_size = local_sites.size();
        neighbor = myrank - 1;
        
        local_sites.resize(org_size+n);
        requests.emplace_back();
        MPI_Irecv(&local_sites[org_size], n, MPI_INT, neighbor, neighbor, MPI_COMM_WORLD, &requests.back());
        a -= n;
        sites_per_rank[myrank] += n;

    }

    //send right
    else if ( b > (target[myrank+1]-1))
    {

        n = std::min((b-target[myrank+1]-1), sites_per_rank[myrank]);
        neighbor = myrank + 1;
        requests.emplace_back();

        org_size = local_sites.size();
        MPI_Isend(&local_sites[0], n, MPI_INT, neighbor, myrank, MPI_COMM_WORLD, &requests.back());
        b -=  n;
        sites_per_rank[myrank] -= n;

    }

    int n_request = requests.size();
    MPI_Waitall(n_request, requests.data(), MPI_STATUSES_IGNORE);


}


