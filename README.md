Here are the notes on K_L identification. It goes through how the variables related to K_L Clusters which were defined in ```DataWriterModule.h``` are defined.

* ```m_KLMnCluster``` - Number of clusters  
  - TYPE: ```Float_t```
  - Defined at ```reconstruction/modules/KlId/KLMExpert/KLMExpertModule.cc```:
    
       ```cpp
       m_KLMnCluster = m_klmClusters.getEntries();
       ```  
    as the number of entries of ```m_klmClusters```.

* ```m_KLMnLayer``` - number of layers hit in KLM cluster
  - TYPE: ```Float_t```
  - Defined at ```reconstruction/modules/KlId/KLMExpert/KLMExpertModule.cc```:
        
      ```cpp
      m_KLMnLayer = cluster.getLayers();
      ```
    as the numbers of layers extracted from function ```getLayers()```, which is defined at ```klm/dataobjects/bklm/BKLMHit1d.h``` like this:
    
      ```cpp
      return BKLMElementNumbers::getLayerByModule(m_ModuleID);
      ```
    where ```BKLMElementNumbers``` is a class defined at
    ```klm/dataobajects/bklm/BKLMElementNumbers.h```.  

  - ```getLayerByModule()``` is also defined there as:
        
      ```cpp
      return ((module & BKLM_LAYER_MASK)>> BKLM_LAYER_BIT) + 1;
      ```

> **Note** CONTINUE LATER!!!!!!!!!!! DONT UNDERSTAND


* ```m_KLMnInnermostLayer```- number of innermost layers hit cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
    
      ```cpp
      m_KLMnInnermostLayer = cluster.getInnermostLayer();
      ```

> **Note** CANT FIND DEFINITION OF getInnermostLayer(), SUSPECT HAS SOMETHING TO DO WITH m\_KLMnLayer, eg the smallest number or sth

    
* ```m_KLMglobalZ```- global Z position in KLM
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as:
    
      ```cpp
      const ROOT::Math::XYZVector& clusterPos = cluster.getClusterPosition();
      m_KLMglobalZ  = clusterPos.Z();
      ```  
    where ```getClusterPosition()``` is defined at ```mdst/dataobjects/ECLClusters.cc```:
    
      ```cpp
      TMatrixDSym ECLCluster::getCovarianceMatrix3x3() const{ 
          const double cluster_x =  getR() * sin(getTheta()) * cos(getPhi())
          const double cluster_y =  getR() * sin(getTheta()) * sin(getPhi());  
          const double cluster_z =  getR() * cos(getTheta());
      return ROOT::Math::XYZVector(cluster_x, cluster_y, cluster_z);
      ```  

> **Note** IS IT WORRYING THAT THE CLUSTER DEFINED HERE IS THE ECLCLUSTER BUT WE ARE USING IT FOR KLM CLUSTER??

> **Note** Cannot find definition of .Z() but I assumed it was extracting the z component of XYZVector

* ```m_KLMtime```- timing of KLM Cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMtime = cluster.getTime();
      ```
> **Note** COULD NOT FIND DEF OF getTime()!!!!

* ```m_KLMavInterClusterDist```- average distance between all KLM clusters
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      tuple<const KLMCluster*, double, double> closestKLMAndDist = findClosestKLMCluster(clusterPos);
      m_KLMavInterClusterDist = get<2>(closestKLMAndDist);
      ```
    where ```findClosestKLMCluster``` is defined at ```reconstruction/modules/KlId/KLMExpert/KlId.h``` (too long to display here).

    From the def of ```findClosestKLMCluster```, which returns:
      ```cpp
      return std::make_tuple(closestKLM, closestKLMDist, avInterClusterDist);
      ```
    it seems reasonable to assume ```get<2>``` extract the 3rd component of ```closestKLMAndDist``` which is the ```m_KLMavInterClusterDist.

    

    









    




