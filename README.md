Here are the notes on K_L identification. It goes through how the variables related to K_L Clusters which were defined in ```DataWriterModule.h``` are defined.

* ```m_KLMnCluster``` - Number of clusters  
  - TYPE: ```Float_t```
  - Defined at ```reconstruction/modules/KlId/KLMExpert/KLMExpertModule.cc```:  
    &nbsp;&nbsp;&nbsp;&nbsp;```m_KLMnCluster = m_klmClusters.getEntries();```  
as the number of entries of ```m_klmClusters```.

* ```m_KLMnLayer``` - number of layers hit in KLM cluster
  - TYPE: ```Float_t```
  - Defined at ```reconstruction/modules/KlId/KLMExpert/KLMExpertModule.cc```:
&nbsp;&nbsp;&nbsp;&nbsp; ```m_KLMnLayer = cluster.getLayers();```     
    
    






