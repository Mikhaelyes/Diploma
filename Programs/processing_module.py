"""Модуль эстиматоров включает в себя:
    -Hill's estimator
    -Ratio estimator
    -Moment estimator
    -UH estimator
    -Pickands estimator
    -Mixed moment

    А так же сторонние функции удобного вывода графиков:
    -plot_est    
    """


import random
import networkx as nx
import pandas as pd
import numpy as np


#Hills estimator
def hill(x):
    """123"""
    maxk = len(x) - 1
    gammah = np.zeros(maxk)
    xord = np.sort(x, axis=0)
    for k in range (0, maxk):
        gammah[k] = (1/(k+1)) * np.sum(np.log(xord[-1:-(k+2):-1])) - np.log(xord[-(k+2)])
    return gammah


#Ratio estimator
def ratio_estimator(x):
    """123"""
    xord = np.sort(x, axis=0)
    level = np.linspace(xord[0], xord[-1], 100)
    rest = np.zeros(len(level))

    for n in range (1, len(level)-1):
        sumup = 0
        sumdn = 0
        for i in x:
            if (i > level[n]):
                sumup += np.log(i / level[n])
                sumdn += 1
        rest[n] = sumup / sumdn
    return rest


#Moment estimator
def moment_estimator(x):
    """123"""
    maxk = len(x) - 1
    xord = np.sort(x, axis=0)
    gammam = np.zeros(maxk)
    gammah = hill(x)
    for k in range (0, maxk):
        sumdn = 0
        for i in range(0, k+1):
            sumdn += (np.log(xord[-1-i+1]) - np.log(xord[-1-k])) ** 2
        sumdn /= k+1
        gammam[k] = gammah[k] + 1 - 0.5 / (1 - (gammah[k] ** 2) / sumdn)
    return gammam


#UH estimator
def uh_estimator(x):
    """123"""
    maxk = len(x) - 1
    gammauh = np.zeros(maxk)
    xord = np.sort(x, axis=0)
    uh = np.zeros(maxk)
    gammah = hill(x)
    for i in range(0, maxk):
        uh[i] = xord[-(i+2)] * gammah[i]
    for k in range (0, maxk-1):
        gammauh[k] = (1/(k+1)) * np.sum(np.log(uh[0:k+1])) - np.log(uh[k+1])
    gammauh[maxk-1] = gammauh[maxk-2]
    return gammauh


#Pickands's estimator
def pickands_estimator(x):
    """123"""
    maxk = len(x) - 1
    xord = np.sort(x, axis=0)
    gammap = np.zeros(round(maxk/4))
    for k in range(0, round(maxk/4)):
        gammap[k] = (1/np.log(2)) * (np.log((xord[-(1+k)] - xord[-(2*k+2)]) \
                                        /(xord[-(2*k+2)] - xord[-(4*k+4)])))
    return gammap


def mixed_moment(x):
    """123"""
    maxk = len(x) - 1
    num_gamma = np.zeros(maxk)
    gammamm_k = 0
    l_mm_n = 0
    m_mm_n = 0
    mm_phi = 0
    xord = np.sort(x, axis=0)
    for k in range(0, maxk):
        l_mm_n = 1 - (1/(k+1)) * np.sum(np.divide(xord[-(k+2)], xord[-1:-(k+2):-1]))
        m_mm_n = (1/(k+1)) * np.sum(np.log(xord[-1:-(k+2):-1])) - np.log(xord[-(k+2)])
        mm_phi = (m_mm_n - l_mm_n) / (l_mm_n * l_mm_n)
        gammamm_k = (mm_phi - 1) / (1 + 2 * np.min((mm_phi - 1), 0))
        num_gamma[k] = gammamm_k
    return num_gamma


def eye_ball(x1):
    """123"""
    x = x1[:int(len(x1)/2)]
    window = int(len(x) / 20)
    diff_x = np.diff(x)
    sum_x = pd.Series()
    for i in range(0, len(x)-window):
        sum_x[i] = np.abs(np.sum(diff_x[i:i+window]))
    k_find = np.argmin(sum_x)
    return x1[k_find], k_find


"""Модуль эстиматоров включает в себя:
    -Hill's estimator
    -Ratio estimator
    -Moment estimator
    -UH estimator
    -Pickands estimator
    -Mixed moment

    А так же сторонние функции удобного вывода графиков:
    -plot_est    
    """


def teil_index_by_sec(df_pr):
    """123"""
    sec = []
    num_iterations = 30
    sec_num = 300
    list_index = pd.Series()
    node_list = list(df_pr['Node'])
    tail_num = 0
#    print(node_list)
    for j in range(sec_num + num_iterations):
        num = random.choice(node_list)
        node_list.remove(num)
#        print(node_list)
        sec.append(df_pr[df_pr['Node'] == num]['PageRank'])
    for i in range(num_iterations):
        tail_data = mixed_moment(sec[i:sec_num+i])
        list_index[i], tail_num = eye_ball(tail_data)
    return list_index.mean(), list_index.std()


def value_index_time(data_vt, comm_list, time):
    """123"""
    # Time_start = int(data_vt['Time_start'].min())
    # Time_stop = int(data_vt['Time_start'].max())
    # step = (Time_stop - Time_start) / 10
    list_time_mean = []
    list_time_std = []
    num_nodes = []
    g1 = nx.Graph()
    for i in time:
        g1.clear()
        df_iter = data_vt[data_vt['Time_start'] < i]
        data_list_t = df_iter[['Source', 'Target']].values.tolist()
        g1 = nx.Graph()
        g1.add_edges_from(data_list_t)
        g1.remove_edges_from(nx.selfloop_edges(g1))
#        print(nx.number_of_nodes(g1))
        pr_t = nx.pagerank(g1, alpha=0.9)
#        print(pr_t)
        df_pr_t = pd.DataFrame(list(pr_t.items()), columns=['Node', 'PageRank'])
#        print(df_pr_t.size)
        df_pr_t_1 = df_pr_t.loc[df_pr_t['Node'].isin(comm_list)]
#        print(df_pr_t_1)
#        df.loc[df['col1']. isin([value1, value2, value3, ...])]
        list_time_mean.append(teil_index_by_sec(df_pr_t_1)[0])
        list_time_std.append(teil_index_by_sec(df_pr_t_1)[1])
        num_nodes.append((df_pr_t_1.size)/2)
    return list_time_mean, list_time_std, num_nodes


def test_tail_index(x):
    """123"""
    global_est = pd.Series()
    list_tail_ind = pd.Series()
    list_step = pd.Series()
    n = len(x)
    h = 300 / n
    step = int(n*h)
    num_blocks = int(1 / h)
    for i in range(num_blocks):
        tail_data = mixed_moment(x[i*step:(i+1)*step])
#        list_tail_ind[i] = Tail_data[50:-50].mean()
        list_tail_ind[i], tail_num = eye_ball(tail_data)

    for i in range(num_blocks):
        global_est[i] = np.sum(list_tail_ind[:i])
        list_step[i] = i/(num_blocks-1)
    return list_step, global_est


def phillips_loretan(x1, x2):
    """123"""
    tail_ind_1, tail_num_1 = eye_ball(mixed_moment(x1))
    tail_ind_2, tail_num_2 = eye_ball(mixed_moment(x2))
    s = (tail_num_1*(tail_ind_2)**2 * (tail_ind_1/tail_ind_2 - 1)**2)/ \
    (tail_ind_1**2 + (tail_num_1/tail_num_2)*tail_ind_2**2)
    return s