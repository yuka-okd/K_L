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

  - ```getLayerByModule()``` is also defined there as:  
```return ((module & BKLM_LAYER_MASK)>> BKLM_LAYER_BIT) + 1;```

CONTINUE LATER!!!!!!!!!!! DONT UNDERSTAND


* ```m_KLMnInnermostLayer```- number of innermost layers hit cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```
```m_KLMnInnermostLayer = cluster.getInnermostLayer();```

CANT FIND DEFINITION OF getInnermostLayer(), SUSPECT HAS SOMETHING TO DO WITH m\_KLMnLayer, eg the smallest number or sth

    
* ```m_KLMglobalZ```- global Z position in KLM
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as:
```const ROOT::Math::XYZVector& clusterPos = cluster.getClusterPosition();```  
```m_KLMglobalZ  = clusterPos.Z();```  
where ```getClusterPosition()``` is defined at ```mdst/dataobjects/ECLClusters.cc```:  
```cpp
TMatrixDSym ECLCluster::getCovarianceMatrix3x3() const{ 
    const double cluster_x =  getR() * sin(getTheta()) * cos(getPhi())
    const double cluster_y =  getR() * sin(getTheta()) * sin(getPhi());  
    const double cluster_z =  getR() * cos(getTheta());
return ROOT::Math::XYZVector(cluster_x, cluster_y, cluster_z); cpp```  

IS IT WORRYING THAT THE CLUSTER DEFINED HERE IS THE ECLCLUSTER BUT WE ARE USING IT FOR KLM CLUSTER??

Cannot find definition of .Z() but I assumed it was extracting the z component of XYZVector





