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
as the numbers of layers extracted from function ```getLayers()```, which is defined at ```klm/dataobjects/bklm/BKLMHit1d.h``` like this:  
```return BKLMElementNumbers::getLayerByModule(m_ModuleID);```
where ```BKLMElementNumbers``` is a class defined at
```klm/dataobajects/bklm/BKLMElementNumbers.h```.  

```getLayerByModule()``` is also defined there as:  
```return ((module & BKLM_LAYER_MASK)>> BKLM_LAYER_BIT) + 1;```

CONTINUE LATER!!!!!!!!!!! DONT UNDERSTAND
  
    
    






