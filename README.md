Here are the notes on K_L identification. It goes through how the variables related to KLM Clusters (defined in ```DataWriterModule.h```) are defined.

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

> **Note** CONTINUE LATER


* ```m_KLMnInnermostLayer```- number of innermost layers hit cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
    
      ```cpp
      m_KLMnInnermostLayer = cluster.getInnermostLayer();
      ```

> **Note** CANT FIND DEFINITION OF getInnermostLayer(), SUSPECT HAS SOMETHING TO DO WITH m_KLMnLayer, eg the smallest number or sth

    
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
> **Note** COULD NOT FIND DEF OF getTime()

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
      double getDistance() const { return m_Distance; }
      ```
    where ```m\_Distance``` is initialized at ```tracking/dataobjects/TrackClusterSeparation.cc```:

      ```cpp
      TrackClusterSeparation::TrackClusterSeparation() :
      m_Distance(1.0E10), // "infinity"
      ```
    > **Note** looks like a dead end but not: ```getSOMETHING``` functions usually grabs the variable that ```setSOMETHING``` has defined.

    Looking at ```tracking/dataobjects/TrackClusterSeparation.h```:
      ```cpp
      void setDistance(double d) { m_Distance = d; }
      ```
    sets the m_Distance that the function ```getDistance()``` retrieves. ```setDistance``` is defeined at ```tracking/trackExtrapolateG4e/TrackExtrapolateG4e.cc``` and the relevant lines are:
      ```cpp
      G4ThreeVector pos  = track->GetPosition(); // this is at postStepPoint
      if (klmClusterInfo != nullptr) {
        for (unsigned int c = 0; c < klmClusterInfo->size(); ++c) {
          G4ThreeVector klmPos = (*klmClusterInfo)[c].second;
          G4ThreeVector separation = klmPos - pos;
          double distance = separation.mag();
          if (distance < klmHit[c].getDistance()) {
            klmHit[c].setDistance(distance);
      ```      
      
> **Note** Cannot find code that changes this value, root files does have varying numbers so it is done somewhere:  ```tracking/trackExtrapolateG4e/TrackExtrapolateG4e.cc``` where? try again

* ```m_KLMTrackSepAngle```- angular distance from track separation object
  - TYPE: ```Float_t```
  - This variable defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      auto trackSeperations = cluster.getRelationsTo<TrackClusterSeparation>();
      trackSep = &trackSeperation;
      m_KLMTrackSepAngle = trackSep->getTrackClusterAngle();
      ```
    where ```getTrackClusterAngle()``` is defined at ```tracking/dataobjects/TrackClusterSeparation.h```:
      ```cpp
      double getTrackClusterAngle() const { return m\_TrackClusterAngle;
      ```
    where again, ```m_TrackClusterAngle``` is initialized in ```tracking/dataobjects/TrackClusterSeparation.cc```:
      ```cpp
      TrackClusterSeparation::TrackClusterSeparation() :
      m_TrackClusterAngle(0.0),
      ```
> **Note** No code was found to alter this value. MUST FIND


* ```m_KLMInitialTrackSepAngle```: angular distance from track to cluster at track starting point
  - TYPE: ```Float_t```
  - The variable is defined at ```reconstruction/modules/KlId/DataWriter/src/DataWriterModule.cc```:
      ```cpp
      m_KLMInitialTrackSepAngle = trackSep->getTrackClusterInitialSeparationAngle();
      ```
    where ```getTrackClusterInitialSeparationAngle()``` is defined at ```tracking/dataobjects/TrackClusterSeparation.h```:

      ```cpp
      double getTrackClusterInitialSeparationAngle() const { return m_TrackClusterInitialSeparationAngle; }
      ```
    Again, initialized at ```tracking/dataobjects/src/TrackClusterSeparation.cc``` but cannot find code that changes this value.
    
    
* ```m_KLMTrackRotationAngle```: angle between track at poca and trackbeginning
  - TYPE: ```Float_t```
  - This varible is defined at  ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMTrackRotationAngle   = -999;
      auto trackSeperations = cluster.getRelationsTo<TrackClusterSeparation>();
      TrackClusterSeparation* trackSep;
      float best_dist = 100000000;
      for (auto trackSeperation :  trackSeperations) {
        float dist = trackSeperation.getDistance();
        if (dist < best_dist) {
          m_KLMTrackRotationAngle = trackSep->getTrackRotationAngle();
      ```
    where ```getTrackRotationAngle()``` is defined at ```tracking/dataobjects/TrackClusterSeparation.h```:
      ```cpp
      double getTrackRotationAngle() const { return m_TrackRotationAngle; }
      ```
    again, ```m_TrackRotationAngle``` is defined at ```tracking/dataobjects/TrackClusterSeparation.cc```:
      ```cpp
      TrackClusterSeparation::TrackClusterSeparation() :
        m_TrackRotationAngle(0.0)
      ```

> **Note** again cannot find code to vary this value, FIND IT
    
      
* ```m_KLMTrackClusterSepAngle```- angle between track momentum and cluster (measured from ip)
  - TYPE: ```Float_t```
  - This variable is defined in a similar way to the previous few variables at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMTrackClusterSepAngle = -999;
      auto trackSeperations = cluster.getRelationsTo<TrackClusterSeparation>();
      for (auto trackSeperation :  trackSeperations) {
        float dist = trackSeperation.getDistance();
        if (dist < best_dist) {
          m_KLMTrackClusterSepAngle = trackSep->getTrackClusterSeparationAngle();
      ```
> **Note** Same as the previous note.

* ```m_KLMAngleToMC```- angle between KLMcluster and Mcparticle
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMAngleToMC = ROOT::Math::VectorUtil::Angle(clusterPos, part->getMomentum());
      ```

> **Note** CANNOT FIND CODE FOR ```Angle()```... Maybe its the same as the function ```ANGLE(P,Q)``` defined at ```generators/koralw/koralw/koeww/kinelib.f```.


* ```m_KLMECLDist```- distance associated ECL <-> KLM cluster
  - TYPE: ```FLoat_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMECLDist = get<1>(closestECLAndDist);
      ```
    where ```closestECLAndDist``` is defined at the same place as :
      ```cpp
      pair<ECLCluster*, double> closestECLAndDist = findClosestECLCluster(clusterPos, eclHypothesis);
      ```
    where ```findClosestECLCluster()``` is defined at ```reconstruction/modules/KlId/KLMExpert/KlId.h``` (long so will not display here).

> **Note** NOW LOOK FOR HOW ```findClosestECLCluster()``` defines ```m_KLMECLDist```


* ```m_KLMECLE```- energy measured in associated ECL cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMECLE = closestECLCluster->getEnergy(eclHypothesis)
      ```
    where ```closestECLCluster``` is, again:
      ```cpp
      pair<ECLCluster*, double> closestECLAndDist = findClosestECLCluster(clusterPos, eclHypothesis);
      ECLCluster* closestECLCluster = get<0>(closestECLAndDist);
      ```
    

> **Note** NOT DONE COME BACK TO THIS




    
    
    

  

    
    
    









    




