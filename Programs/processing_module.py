"""
Модуль обработки предназначен для проведения исследования последовательностей и графов.
Состоит из следующих частей:
    1) Модуль эстиматоров.
    2) Модуль оценки стационарности последовательностей.
    3) Модуль оценки стационарности графов.
"""


import random
import networkx as nx
import pandas as pd
import numpy as np
from scipy import stats


"""

МОДУЛЬ ЭСТИМАТОРОВ 

Включает в себя:
    -Hill's estimator
    -Ratio estimator
    -Moment estimator
    -UH estimator
    -Pickands estimator
    -Mixed moment

    А так же функции:
    -eye_ball - оценки плато.
    -bootstrap_est - рассчёта доверительного интервала для эстиматоров.
"""


def hill(x):
    """
    Функция позволяет сделать оценку Хилла индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        gammah - оценка.
    """

    maxk = len(x) - 1
    gammah = np.zeros(maxk)
    xord = np.sort(x, axis=0)
    for k in range (0, maxk):
        gammah[k] = (1/(k+1)) * np.sum(np.log(xord[-1:-(k+2):-1])) - np.log(xord[-(k+2)])
    return gammah


def ratio_estimator(x):
    """
    Функция позволяет сделать оценку Ratio индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        rest - оценка.
    """

    xord = np.sort(x, axis=0)
    level = np.linspace(xord[0], xord[-1], 100)
    rest = np.zeros(len(level))

    for n in range (1, len(level)-1):
        sumup = 0
        sumdn = 0
        for i in x:
            if(i > level[n]):
                sumup += np.log(i / level[n])
                sumdn += 1
        rest[n] = sumup / sumdn
    return rest


def moment_estimator(x):
    """
    Функция позволяет сделать оценку moment_estimator индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        gammam - оценка.
    """

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


def uh_estimator(x):
    """
    Функция позволяет сделать оценку uh_estimator индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        gammauh - оценка.
    """

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


def pickands_estimator(x):
    """
    Функция позволяет сделать оценку pickands_estimator индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        gammap - оценка.
    """

    maxk = len(x) - 1
    xord = np.sort(x, axis=0)
    gammap = np.zeros(round(maxk/4))
    for k in range(0, round(maxk/4)):
        gammap[k] = (1/np.log(2)) * (np.log((xord[-(1+k)] - xord[-(2*k+2)]) \
                                        /(xord[-(2*k+2)] - xord[-(4*k+4)])))
    return gammap


def mixed_moment(x):
    """
    Функция позволяет сделать оценку mixed_moment индекса экстримальной величины.
    Input:
        x - последовательность.
    Output:
        num_gamma - оценка.
    """

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


def eye_ball(x):
    """
    Функция реализует метод eye ball для некоторой функции.
    Input:
        x - последовательность определённая некоторой функцией.
    Output:
        out - Значение оценщика по методу eye ball.
        k_find + int(w/2) - значение порядковой статистики соответствующей значению оценщика.
    """

    w = int(len(x)/50)
    n = len(x)
    e = 0
    h = 0.9
    list_k = []
    eye_int = 0
    e = 0.1 * np.mean(x)
    # print(e)
    # for i in range(int(0.15*n), n-w):
    for k in range(2, n-w):
        # e = 0.1 * x[i]
        eye_int = 0
        for i in range(w):
            if (np.abs(x[k+i] - x[k]) < e):
                y = 1
            else:
                y = 0
            eye_int += y
        if (eye_int/w > h):
            list_k.append(k)
    k_find = np.min(list_k)
    out = x[k_find + int(w/2)]
    return out, k_find + int(w/2)


def bootstrap_est(x, interval=0.9, func=None):
    """
    Функция позволяет вычислить доверительный интервал для некоторого эстиматора.
    Input:
        x - последовательность.
        interval - квантиль доверительного интервала.
        func - применяемый эстиматор.
    Output:
        output_mn - среднее.
        output_up - верхняя грань доверительного интервала.
        output_dn - нижняя грань доверительного интервала.
    """

    len_x = len(x)
    data = list(x)
    num_b = 100

    y = np.zeros((num_b, len_x))
    gamma_est = np.zeros((num_b, len(func(x))))    
    for i in range(num_b):
        for j in range(len_x):
            y[i,j] = random.choice(data)
        gamma_est[i,:] = func(y[i,:])

    len_seq = len(gamma_est[1,:])
    output = np.zeros((3, len_seq))
    for i in range(len_seq-1):
        interval_e = stats.norm.interval(interval, loc=np.mean(gamma_est[:, i]), scale=gamma_est[:, i].std())
        output[0, i] = gamma_est[:, i].mean()
        output[1, i] = interval_e[0]
        output[2, i] = interval_e[1]

    return output


"""

МОДУЛЬ ГЕНЕРАЦИИ ГРАФОВ 

Включает в себя:
    -gen_graph_PA
    -gen_graph_CA
    -gen_graph_ABG
"""


def gen_graph_PA(G_input: nx.Graph, num_nodes: int, num_neigh: int) -> nx.Graph:
    """
    Функция генерации вершин графа методом предпочтительного присоединения.
    Используются библиотеки tqdm, networkx.
    
    Input:
        G_input - изначальный неориентированный граф к которому присоединяются новые вершины.
        num_nodes - количество присоединяемых вершин.
        num_neigh - количество соседий, к которым присоединяется новая вершина.
        
    Output:
        Изначальный неориентированный граф к которому присоеденены вершины.
    """
    i_start = len(list(G_input.nodes))
    print(i_start)
    for i in tqdm(range(i_start + 1, i_start + num_nodes)):
        k = []
        for i1 in range(num_neigh):
            k.append((random.choices(list(dict(G_input.degree()).keys()), weights=list(dict(G_input.degree()).values())))[0])
        for i1 in range(num_neigh):
            G_input.add_edge(k[i1], i)     
    return G_input


def gen_graph_CA(G_input: nx.Graph, num_nodes: int, num_neigh: int) -> nx.Graph:
    """
    Функция генерации вершин графа методом кластерного присоединения.
    Используются библиотеки tqdm, networkx.
    Содержит циклы так как k.append((random.choices)) может выдавать 2 одинаковых числа. 
    
    Input:
        G_input - изначальный неориентированный граф к которому присоединяются новые вершины.
        num_nodes - количество присоединяемых вершин.
        num_neigh - количество соседий, к которым присоединяется новая вершина.
        
    Output:
        Изначальный неориентированный граф к которому присоеденены вершины.
    """
    i_start = len(list(G_input.nodes))
    for i in tqdm(range(i_start + 1, i_start + num_nodes)):
        local_calst = list(nx.clustering(G_input).values())
        for j in range(len(local_calst)):
            local_calst[j] += 0.001
        k = []
        for i1 in range(num_neigh):
            k.append((random.choices(list(nx.clustering(G_input).keys()), weights=local_calst))[0])
        for i1 in range(num_neigh):
            G_input.add_edge(k[i1], i)
    return G_input


def gen_graph_ABG(G_input: nx.DiGraph, num_iterations: int, alpha: float, beta: float, \
                  d_in: float, d_out: float) -> nx.DiGraph:
    """
    Функция генерации вершин графа методом альфа, бетта, гамма присоединения.
    Используются библиотеки tqdm, networkx.
    gamma = 1 - (alpha + beta).
    
    Input:
        G_input - изначальный ориентированный граф к которому присоединяются новые вершины.
        num_iterations - количество итераций алгоритма.
        alpha - float вероятность.
        beta - float вероятность. 
        d_in - float.
        d_out - float.
        
    Output:
        Изначальный граф к которому присоеденены вершины.
    """
    assert ((alpha + beta) <= 1)

    i_start = len(list(G_input.nodes))
    gamma = 1 - alpha - beta
    print(alpha, beta, gamma, G_input)

    for i in tqdm(range(0, num_iterations + 1)):
        iter_prob = random.choices([1, 2, 3], weights=[alpha, beta, gamma])[0]
        
        I_n1 = list(dict(G_input.in_degree()).values())
        O_n1 = list(dict(G_input.out_degree()).values())
        
        N_I = list(dict(G_input.in_degree()).keys())
        N_O = list(dict(G_input.out_degree()).keys())
        
        N = len(I_n1)
        n = sum(list(dict(G_input.in_degree()).values()))
        u = list(G_input.nodes)[-1] + 1

        assert (N_I == N_O)
        
        if (iter_prob == 1):
            denominator = n - 1 + d_in * N
            local_eq = []
            for j in range(len(I_n1)):
                local_eq.append((I_n1[j] + d_in) / denominator)
            w = (random.choices(N_I, weights=local_eq))[0]
            G_input.add_edge(i, w)

        if (iter_prob == 2):
            denominator_1 = n - 1 + d_in * N
            denominator_2 = n - 1 + d_out * N
            local_eq = []
            for j in range(len(I_n1)):
                local_eq.append(((I_n1[j] + d_in) / denominator_1) * ((O_n1[j] + d_out) / denominator_2))
            w = (random.choices(N_I, weights=local_eq))[0]
            u = (random.choices(N_I, weights=local_eq))[0]
            G_input.add_edge(u, w)
        
        if (iter_prob == 3):
            denominator = n - 1 + d_out * N
            local_eq = []
            for j in range(len(O_n1)):
                local_eq.append((O_n1[j] + d_out) / denominator)
            w = (random.choices(N_O, weights=local_eq))[0]
            G_input.add_edge(w, i)

    return G_input


"""

МОДЕЛЬ ОЦЕНКИ СТАЦИОНАРНОСТИ ПОСЛЕДОВАТЕЛЬНОСТЕЙ

Включает в себя:
    -test_tail_index
    -phillips_loretan
"""


def test_tail_index(x, amount=300):
    """
    Функция позволяет найти глобальный оценщик.
    """
    global_est = pd.Series()
    list_tail_ind = pd.Series()
    list_step = pd.Series()
    n = len(x)
    h = amount / n
    step = int(n*h)
    num_blocks = int(1 / h)
    for i in range(num_blocks-1):
        tail_data = hill(x[i*step:(i+1)*step])
#        list_tail_ind[i] = Tail_data[50:-50].mean()
        list_tail_ind[i], tail_num = eye_ball(tail_data)

    for i in range(num_blocks):
        global_est[i] = np.sum(list_tail_ind[:i])
        list_step[i] = i/(num_blocks-1)
    return list_step, global_est


def phillips_loretan(x1, x2):
    """
    Функция теста phillips_loretan для двух последовательностей
    """
    tail_ind_1, tail_num_1 = eye_ball(mixed_moment(x1))
    tail_ind_2, tail_num_2 = eye_ball(mixed_moment(x2))
    s = (tail_num_1*(tail_ind_2)**2 * (tail_ind_1/tail_ind_2 - 1)**2)/ \
    (tail_ind_1**2 + (tail_num_1/tail_num_2)*tail_ind_2**2)
    return s


"""

МОДУЛЬ ОЦЕНКИ СТАЦИОНАРНОСТИ ГРАФОВ

Включает в себя:
    -teil_index_by_sec
    -value_index_time
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