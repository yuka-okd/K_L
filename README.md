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
    it seems reasonable to assume ```get<2>``` extract the 3rd component of ```closestKLMAndDist``` which is the ```m_KLMavInterClusterDist```.


* ```m_KLMhitDepth```- hit depth in KLM, distance to IP
  - TYPE: ```Float_t```
  - This variable defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMhitDepth = cluster.getClusterPosition().R();
      ```
> **Note** COOULD NOT FIND THE DEFINITION OF R(), PERHAPS FIND RADIUS?

* ```m_KLMenergy```- Energy deposit in KLM (0.2 GeV * nHitCells)
  - TYPE: ```Float_t```
  - This variable defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m\_KLMenergy  = cluster.getEnergy();
      ```
    where ```getEnergy()``` is defined at: ```analysis/dataobjects/Particle.h```:
      ```cpp
      return sqrt(m_momentumScale * m_momentumScale * m_px * m_px + m_momentumScale * m_momentumScale * m_py * m_py +  m_momentumScale * m_momentumScale * m_pz * m_pz + m_mass * m_mass) - m_energyLossCorrection;

> **Note** Continue: find the variables in the squreroot equation

* ```m_KLMinvM```- invariant mass calculated from root vector
  - TYPE: ```Float_t```
  - This variable defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMinvM  = cluster.getMomentum().M2();
      ```
    where ```getMomentum``` is defined at ```mdst/dataobjects/MCParticle.h```:
      ```cpp
      ROOT::Math::XYZVector getMomentum() const { return ROOT::Math::XYZVector(m_momentum_x, m_momentum_y,
      m_momentum_z); }
      ```

* ```m_KLMTruth```- target variable for KLM classification
  - TYPE: ```Float_t```
  - defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMTruth = mcParticleIsKlong(part);
      ```
    where ```mcParticleIsKlong``` is defined at ```reconstrucion/KlId/KLMExpert/KlId.h``` (code
    is long, will not write here) this function will return ```true``` if the input particle is K_L
    or if the particle is a daughter of K_L. Otherwise will return ```false```. 

> **Note** this function returns 0 or 1. Is it reasonable to define ```m\_KLMTruth``` as ```Float_t``` rather than ```int``` type??

    
* ```m_KLMnextCluster```- distance to next KLM cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMnextCluster = get<1>(closestKLMAndDist);
      ```
      As mentioned before, ```closestKLMAndDist``` is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as the output of function ```findClosestKLMCluster```, which is defined at ```recontruction/modules/KlId/KLMExpert/KlId.h```:
      ```cpp
      std::make_tuple(closestKLM, closestKLMDist, avInterClusterDist)
      ```
    Therefore, ```m_KLMnextCluster``` is defined as the 2nd component of ```closestKLMAndDist```, which is ```closestKLMDist```.

* ```m_KLMTrackSepDist```- distance from track separation object
  - TYPE: ```Float_t```
  - This variable is defined at: ```reconstruction/modules/KlId/KLMExpert/KLMExpertModule.cc```
      ```cpp
      m_KLMTrackSepDist = -999;
      auto trackSeperations = cluster.getRelationsTo<TrackClusterSeparation>();
      if (dist < best_dist) {
      best_dist = dist
      m_KLMTrackSepDist = trackSeperation.getDistance();
      ```
    where ```getRelationsTo()``` is defined at ```framework/datastore/RelationsObject.h```:

      ```cpp
      return RelationVector<TO>(DataStore::Instance().getRelationsWith(DataStore::c_ToSide, this, m_cacheDataStoreEntry, m_cacheArrayIndex, TO::Class(), name, namedRelation));
      ```
    and ```getDistance()``` defined at ```tracking/dataobjects/TrackClusterSeparation.h```:

      ```cpp
      double getDistance() const { return m\_Distance; }
      ```
    where ```m\_Distance``` is initialized at ```tracking/dataobjects/TrackClusterSeparation.cc```:

      ```cpp
      TrackClusterSeparation::TrackClusterSeparation() :
      m_Distance(1.0E10), // "infinity"
      ```
> **Note** Cannot find code that changes this value, would be bad if this was used in all cases, FIND WHERE


    
    
      
    
    
    
    
    
    

  

    
    
    









    




