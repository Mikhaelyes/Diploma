## Магистерский диплом ##

1) Data_initial - папка с изначальными данными графов. <br/>
   Содержит: fb-messages.edges. <br/>
2) Data_postproc - папка с изменением данных для определённых целей. <br/>
3) Programs - папка с программами.<br/>
   Содержит:<br/>
       -processing_module - модуль робработки последовательностей.<br/>
       -Real_graph - программа визуализации экспериментов с реальным графом.<br/>
       -Synthetic_sequence - программа визуализации экспериментов с синтетической последовательностью.<br/>
       -Synthetic_graph - программа визуализации экспериментов с случайным графом.<br/>
       -Old_versions - старые версии. <br/>
4) Text - текст магистерского диплома. <br/>

processing_module содержит: <br/>
1) МОДУЛЬ ЭСТИМАТОРОВ <br/>
Включает в себя: <br/>
    -Hill's estimator <br/>
    -Ratio estimator <br/>
    -Moment estimator <br/>
    -UH estimator <br/>
    -Pickands estimator <br/>
    -Mixed moment <br/>
А так же функции: <br/>
    -eye_ball - оценки плато. <br/>
    -bootstrap_est - рассчёта доверительного интервала для эстиматоров. <br/>
2) МОДУЛЬ ГЕНЕРАЦИИ ГРАФОВ  <br/>
Включает в себя: <br/>
    -gen_graph_PA <br/>
    -gen_graph_CA <br/>
    -gen_graph_ABG <br/>
3) МОДЕЛЬ ОЦЕНКИ СТАЦИОНАРНОСТИ ПОСЛЕДОВАТЕЛЬНОСТЕЙ <br/>
Включает в себя: <br/>
    -test_tail_index <br/>
    -phillips_loretan <br/>
4) МОДУЛЬ ОЦЕНКИ СТАЦИОНАРНОСТИ ГРАФОВ <br/>
Включает в себя: <br/>
    -teil_index_by_sec <br/>
    -value_index_time <br/>
