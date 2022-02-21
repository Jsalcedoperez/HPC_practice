#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <vector>

namespace defi{

    int fsc_size, new_size;
    std::vector<int> sites;
    int ierr, myrank, tot_ranks;
    int local_fsc_size;
    int a, b, n, neighbor;
    std::vector<int> min_work(5,0);
    int total_sites, remainder, min_target;
    std::vector<MPI_Request> requests;

}


void setup(const std::vector<int>&);

void create_sites();

void rebalance();

int main(int argc, char *argv[])

{

    using namespace defi;

    std::vector<int> fission_sites = {0, 5, 10, 18, 25};

    setup(fission_sites);

    create_sites();
    
    int ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &defi::myrank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &defi::tot_ranks);


    int tot_converged, converge = 0;

    while (tot_converged < tot_ranks)
    {

        rebalance();
        fsc_size = b - a;

        if (fsc_size == min_work[myrank] )
        {

            converge = 1;

        }

         MPI_Allreduce(&converge, &tot_converged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    }


    ierr = MPI_Finalize();

    return 0;
}

void setup(std::vector<int>& fission_sites)
{
    using namespace defi;
    fsc_size = fission_sites.size();

    local_fsc_size = fission_sites[myrank + 1];
    a = fission_sites[myrank];
    b = a + local_fsc_size - 1;

    total_sites = fission_sites[fsc_size];
    remainder = total_sites % tot_ranks;
    min_target = total_sites / tot_ranks;

    for (int irank = 0 ; irank < tot_ranks ; ++irank)
    {

        int imin_work = irank < remainder ? min_target + 1 : min_target;
        min_work.push_back(imin_work);

    }

}

void create_sites()
{

    //set up sites
    using namespace defi;
    if (sites.empty())
    {

        for (int i = a; i < b; ++i)
        {
            sites.push_back(i);

        }
    }

    else

    {
        new_size = b-a;
        sites.resize(new_size);
        if (new_size > fsc_size){
            for (int i = b ; i < new_size ; ++i)
            {
                sites.push_back(i);
            
        }
    }

    fsc_size = b - a;

    if (fsc_size != sites.size())
    {

        std::cout << "sites' size and fsc_size don't match" << std::endl;

    }

}

void rebalance()

{

    using namespace defi;
    if (min_work[myrank] > a)
    {

        n = std::min((min_work[myrank] - a), local_fsc_size);
        neighbor = myrank - 1;
        requests.emplace_back();
        MPI_Isend(&sites, n, MPI_INT, neighbor, myrank, MPI_COMM_WORLD, &requests.back());
        a += n;

    }

    else if (min_work[myrank + 1] < b)
    {
        n = std::min((b - min_work[myrank+1]), local_fsc_size);
        neighbor = myrank + 1;
        requests.emplace_back();
        MPI_Isend(&sites, n, MPI_INT, neighbor, myrank, MPI_COMM_WORLD, &requests.back());
        b -=  n;


    }

//receive
    if (min_work[myrank + 1] > b)
    {
        n = std::min((min_work[myrank + 1] - b), local_fsc_size);
        neighbor = myrank + 1;
        sites.resize(sites.size() + n - 1);
        requests.emplace_back();
        MPI_Irecv(&sites, n, MPI_INT, neighbor, neighbor, MPI_COMM_WORLD, &requests.back());
        b += n;

    }

    else if (min_work[myrank] < a)
    {

        n = std::min((a - min_work[myrank]), local_fsc_size);
        neighbor = myrank - 1;
        sites.resize(sites.size() + n - 1);
        requests.emplace_back();
        MPI_Irecv(&sites, n, MPI_INT, neighbor, neighbor, MPI_COMM_WORLD, &requests.back());
        a -= n;


    }

    int n_request = requests.size();
    MPI_Waitall(n_request, requests.data(), MPI_STATUSES_IGNORE);


}


