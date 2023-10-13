### Based on https://mikecvet.wordpress.com/2010/07/02/parallel-mapreduce-in-python/
import json
import sys
import igraph
import numpy as np
from multiprocessing import Pool 
import time 
import os
import logging
import datetime
import warnings
import pandas as pd 

sys.path.insert(0, '/Users/bz247')
from sp import interface 

absolute_path = os.path.dirname(os.path.abspath(__file__))

def map_edge_pop(row):
    ### Find shortest path for each unique origin --> one destination
    ### In the future change to multiple destinations

    origin_ID = int(OD['O'].iloc[row]) + 1 ### origin's ID on graph nodes
    destin_ID = int(OD['D'].iloc[row]) + 1 ### destination's ID on graph nodes
    traffic_flow = int(OD['flow'].iloc[row]) ### number of travellers with this OD

    sp = g.dijkstra(origin_ID)
    if sp.distance(destin_ID) > 10e7:
        return [], 0
    sp_route = sp.route(destin_ID)
    results = []
    ### how to do the in-place updating?
    # for edge in sp_route:
    #     graph_v1 = edge[0]-1
    #     graph_v2 = edge[1]-1
    #     results.append([graph_v1, graph_v2, traffic_flow])
    return results, 1

def edge_tot_pop(L, day, hour):
    logger = logging.getLogger('main.one_step.edge_tot_pop')
    t0 = time.time()
    edge_volume = {}
    for sublist in L:
        for p in sublist: ### p[0] is edge ID on graph, p[1] is edge flow
            try:
                edge_volume[p[0]] += p[1]
            except KeyError:
                edge_volume[p[0]] = p[1]
    t1 = time.time()
    logger.info(
        f'DY{day}_HR{hour}: # edges to be updated {len(edge_volume)}, taking {t1 - t0} seconds'
    )

    # with open('edge_volume_2.json', 'w') as outfile:
    #     json.dump(edge_volume, outfile, indent=2)

    return edge_volume

def one_step(day, hour):
    ### One time step of ABM simulation

    logger = logging.getLogger('main.one_step')

    ### Read the OD table of this time step
    global OD
    OD = pd.read_csv(
        absolute_path
        + f'/../TNC/output/SF_graph_DY{day}_HR{hour}_OD_50000.csv'
    )

    ### Number of processes
    process_count = 4
    logger.debug(f'number of process is {process_count}')

    ### Build a pool
    pool = Pool(processes=process_count)
    logger.debug('pool initialized')

    ### Find shortest pathes
    unique_origin = 200 # OD.shape[0]
    logger.info(f'DY{day}_HR{hour}: # OD rows (unique origins) {unique_origin}')
    t_odsp_0 = time.time()
    res = pool.imap_unordered(map_edge_pop, range(unique_origin))

    ### Close the pool
    pool.close()
    pool.join()
    t_odsp_1 = time.time()
    logger.debug(f'shortest_path running time is {t_odsp_1 - t_odsp_0}')

    ### Collapse into edge total population dictionary
    edge_pop_tuples, destination_counts = zip(*res)
    logger.info(f'DY{day}_HR{hour}: # destinations {sum(destination_counts)}')
    #edge_volume = edge_tot_pop(edge_pop_tuples, day, hour)
    #print(list(edge_volume.items())[0])

    #return edge_volume
    return 0

def main():
    logging.basicConfig(filename=absolute_path+'/sf_abm_mp.log', level=logging.DEBUG)
    logger = logging.getLogger('main')
    logger.info(f'{datetime.datetime.now()} \n\n')

    t_start = time.time()

    ### Read initial graph
    #global g_igraph
    #g_igraph = igraph.Graph.Read_Pickle(absolute_path+'/../data_repo/data/sf/network_graph.pkl')

    global g
    g = interface.readgraph(bytes(absolute_path+'/../data_repo/data/sf/network_sparse.mtx', encoding='utf-8'))

    for day in [1]:
        for hour in range(9, 10):

            logger.info(f'*************** DY{day} HR{hour} ***************')

            t0 = time.time()
            edge_volume = one_step(day, hour)
            t1 = time.time()
            logger.info(f'DY{day}_HR{hour}: running time {t1 - t0}')

                    ### Update graph
                    # volume_array = np.zeros(g.ecount())
                    # volume_array[list(edge_volume.keys())] = np.array(list(edge_volume.values()))*400 ### 400 is the factor to scale Uber/Lyft trip # to total car trip # in SF.
                    # g.es['volume'] = volume_array
                    # logger.info('DY{}_HR{}: max link volume {}'.format(day, hour, max(volume_array)))
                    # g.es['t_new'] = fft_array*(1.2+0.78*(volume_array/capacity_array)**4) ### BPR and (colak, 2015)
                    # g.es['weight'] = g.es['t_new']

                    #write_geojson(g, day, hour)

    t_end = time.time()
    logger.info(f'total run time is {t_end - t_start} seconds \n\n\n\n\n')

if __name__ == '__main__':
    main()

