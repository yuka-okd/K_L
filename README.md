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
    where ```findClosestECLCluster``` returns: ```return std::make_pair(closestECL, closestECLAngleDist);```.
    
    

> **Note** UNDERSTAND WHAT `findClosestCluster``` DOES: IT IS DEFINED AT ```KlId.h```

* ```m_KLMECLdeltaL```- distance between track entry point and cluster center, might be removed
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc```:
      ```cpp
      m_KLMECLdeltaL         = closestECLCluster->getDeltaL();;
      ```
    where ```closestECLCluster``` is defined in the same code as
      ```cpp
      pair<ECLCluster*, double> closestECLAndDist = findClosestECLCluster(clusterPos, eclHypothesis);
      ECLCluster* closestECLCluster = get<0>(closestECLAndDist);
      ```
    ```getDeltaL()``` is defined to get the value of ```m_deltaL``` in     
    ```mdst/dataobjects/ECLCluster.h```, and ```m_deltaL``` is defined in ```ecl/modules/eclClusterProperties/ECLClusterPropertiesModules.cc``` like below:

      ```cpp
      // compute path lengths on the energy weighted average crystals
      // direction and on the extrapolated track direction corresponding to
      // the minimum distance among the two lines. if more than one track is
      // related to a cluster the one with the highest momentum is used
      if (cluster->isTrack()) {
        double lTrk, lShower;
        computeDepth(shower, lTrk, lShower);
        cluster->setdeltaL(lTrk);
      ```
    where ```lTrk``` is given a value in funcion ```computeDepth()```.

    ```computeDepth()``` is defined in the same file with relevant snippit as follows:
      ```cpp
      ROOT::Math::XYZVector cvec = geometry->GetCrystalVec(cellid - 1);
      avgDir += energy * cvec; 
      ROOT::Math::XYZVector w0 = showerCenter - trkpos;
      double costh = avgDir.Unit().Dot(trkdir);
      double sin2th = 1 - costh * costh;
  
      lTrk = w0.Dot(trkdir) - costh * w0.Dot(avgDir.Unit());
      lTrk /= sin2th;
      ```
    ```avgDir``` seems to represent the general direction of the shower. ```w0``` is a vector between track position and center of shower.  ```trkdir``` is defined in the same code as ```trkdir = extHit.getMomentum().Unit();``` which seems to be the unit vector direction of the hit. Therefore costh is the dot product between the shower and the track which gives cosine angle between the two. The last two lines calculates the distance between track and shower with respect to the track's direction.


>**Note** WHAT IS THE DIVISION FOR


* ```m_KLMECLminTrackDist```- track distance between associated ECL cluster and track extrapolated into ECL
  - TYPE: ```Float_t```
  - This varible is defined in ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as
      ```cpp
      m_KLMECLminTrackDist   = closestECLCluster->getMinTrkDistance();
      ```
    where the function ```getMinTrkDistance()``` is defined to extract the output of ```setMinTrkDistance(dist)``` which is defined in ```ecl/modules/eclClusterProperties/ECLClusterPropertiesModules.cc``` as ```double dist = computeTrkMinDistance(shower, m_tracks, trackID);``` where ```computeTrkMinDistance()``` is defined below as:
      ```cpp
      ROOT::Math::XYZVector cryCenter;
      VectorUtil::setMagThetaPhi(
      cryCenter, shower.getR(), shower.getTheta(), shower.getPhi());
      for (const auto& extHit : track.getRelationsTo<ExtHit>()) {
        trkpos = extHit.getPosition();
        double distance = (cryCenter - trkpos).R();
      }
      ```
    ```setMaghetaPhi()``` is defined at ```framework/geometry/VectorUtil.h``` where it converts the mag, theta and phi coordinates to x, y and z coordinates and outputs the vector in cartesian. ```cryCenter``` refers to the geometrical center of the shower in ECL. ```trkpos``` extracts the position of the extrapoalted hit in the ECL. ```distance``` equals to the distance between the centre of the shower in ECL and the extrapolated hit in the ECL. This becomes the input value for ```setMinTrkDistance``` which outputs ```m_KLMECLminTrackDist = distance```. In my understanding, the smaller the ```m_KLMECLminTrackDist```, the more accurate the prediction is. 
      
    
 
* ```m_KLMECLE9oE25```- E in surrounding 9 crystals divided by surrounding 25 crydtalls
  - TYPE: ```Float_t```
  - This varibale is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```closestECLCluster->getE9oE21();```, where the output of ```getE90E21()``` is set by the return value of ```setE9oE21()```
    ```setE9oE25()``` is defined at ```ecl/dataobjects/ECLShower``` as ```void setE9oE21(double E9oE21) { m_E9oE21 = E9oE21; }```. The input for the ```setE9oE21``` appears at ```ecl/modules/eclShowerShape/ECLShowerShapeModule.cc`` as shown: ```if (eclShower->getE9oE21() < 1e-9) eclShower->setE9oE21(computeE9oE21(*eclShower));```. ```m_E9oE21``` is initiated at ```ecl/dataobjects/ECLShower.h``` as 0.0. ```computeE9oE21()``` is defined in ```ecl/modules/eclShowerShape/ECLShowerShapeModule.cc```:
      ```cpp
      const auto it9 = std::find(n9.begin(), n9.end(), cellid);
      if (it9 != n9.end()) {
        energy9 += weight * energy;
      }
      const auto it21 = std::find(n21.begin(), n21.end(), cellid);
      if (it21 != n21.end()) {
        energy21 += weight * energy;
      }
      ```
    where ```n9``` and ```n21``` are the list of 9 and 21 neighbour ids. Therefore the function returns the sum of weighted energy of the neighbours. 



>**Note** It says E25 but actually uses E21 which is the 5x5 crystals excluding the four corners 

 
* ```m_KLMECLTiming```- timing of associated ECL cluster
  - TYPE: ```Float_t```
  - This varibale is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMECLTiming         = closestECLCluster->getTime();```, where ```getTime()``` returns ```m_time```.
>**Note** Check whether setTime() is defined at ```rawdata/modules/SetMetaTimeModule.py``` or not 

* ```m_KLMECLTerror```- uncertainty on time in associated ECL cluster
  - TYPE: ```FLoat_t```
  - This variable defined in the ```DataWriterModule.cc``` as ```m_KLMECLTerror         = closestECLCluster->getDeltaTime99();``` where ```getDeltaTime99()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ```double getDeltaTime99() const {return m_deltaTime99;}```, where ```m_deltaTime99``` is is set by ```setDeltaTime99()```, which is defined as tim ethat contains 99% of signal crystals.
    ```m_deltaTime99``` seemed to be defined at ```ecl/modules/eclSplitterN1/ECLSplitterN1Module.cc```:
      ```cpp
      for (unsigned int i = 0; i < digitVector.size(); ++i) {
        const ECLCalDigit dig = digitVector[i];
        highestEnergyTimeResolution = dig.getTimeResolution();
        aECLShower->setDeltaTime99(highestEnergyTimeResolution);
      ```
    where ```highestEnergyTimeResolution``` is defined in ```ecl/modules/eclDigitCalibration/ECLDigitCalibratorModule.cc```:
      ```cpp
      for (auto& aECLCalDigit : m_eclCalDigits) {

        // perform the time resolution calibration
        const double t99 = getT99(aECLCalDigit.getCellId(),
                                  aECLCalDigit.getEnergy(),
                                  aECLCalDigit.hasStatus(ECLCalDigit::c_IsFailedFit),
                                  bgCount); // calibrated time resolution
        aECLCalDigit.setTimeResolution(t99);
      ```
    where ```getT99()``` is defined in the same file as:
      ```cpp
      double ECLDigitCalibratorModule::getT99(const int cellid, const double energy, const bool fitfailed, const int bgcount) const
      {
      const double bglevel = TMath::Min((double)bgcount / (double)c_nominalBG * m_th1fBackground->GetBinContent(cellid) / m_averageBG, m_pol2Max); // c_nominalBG = 183 for actual version of digitizer, m_averageBG is about 2 MeV/ mus
      // Get p1 as function of background level
      const double p1 = c_pol2Var1 + c_pol2Var2 * bglevel + c_pol2Var3 * bglevel * bglevel;

      // inverse energy in 1/MeV
      double einv = 0.;
      if (energy > 1e-12)
        einv = 1. / (energy / Belle2::Unit::MeV);

      if (energy > 1e-12)
        einv = 1. / (energy / Belle2::Unit::MeV);
        // calculate t99 using the inverse energy and p1 (p0 is zero)
        double t99 = p1 * einv;
    
      // for high energies we fix t99 to 3.5ns
      if (t99 < c_minT99)
        t99 = c_minT99;
      ```
<**Note** investigate what p1 and other variables in bglevel is

    
  
* ```m_KLMECLEerror```- uncertainty on E in associated ECL cluster
  - TYPE: ```Float_t```
  - Variable defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMECLEerror         = closestECLCluster->getUncertaintyEnergy();```. ```mdst/dataobjects/ECLCluster.h``` defines ```double getUncertaintyEnergy() const {return (m_sqrtcovmat_00);}```, where ```m_sqrtcovmat_00``` is set at ```mdst/dataobjects/ECLCluster.h``` as ```void setCovarianceMatrix(double covArray[6])
    {
      m_sqrtcovmat_00 = sqrt(fabs(covArray[0]));```. The input matrix is defined at ```ecl/eclCovarianceMatrix/ECLCovariantMatrixModule.cc``` as: ``` double covMatrix[6] = {sigmaEnergy * sigmaEnergy, 0.0, sigmaPhi * sigmaPhi, 0.0, 0.0, sigmaTheta * sigmaTheta};``` where the relevant ```covMatrix[0]``` is in terms of ```sigmaEnergy``` which is defined differently according to the detector region, and ```const double energy = eclShower.getEnergy();``` This script only does the calculation for photon showers. 
    

* ```m_KLMtrackToECL```- primitive distance cluster <-> track for associated ECL cluster
  - TYPE: ```Float_t```
 
<**Note** CANNOT FIND ANY INFORMATION ON THIS

* ```m_KLMKLid```- KlId for that object
  - TYPE: ```Float_t```
  - This varibale is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as:
      ```cpp
      KlId* klid = cluster.getRelatedTo<KlId>();
      if (klid) {
      m_KLMKLid = klid->getKlId();
      ```
    where ```getKlId()``` is defined at ```mdst/dataobjects/KlId.cc```:
      ```cpp
      double KlId::getKlId() const
      {
        auto klmClusterWeight = getRelatedFromWithWeight<KLMCluster>();
        if (klmClusterWeight.first) return klmClusterWeight.second;
        auto eclClusterWeight = getRelatedFromWithWeight<ECLCluster>();
        if (eclClusterWeight.first) return eclClusterWeight.second;
        return nan("");
      }
      ```
    where ```getRelatedFromWithWeight()``` is defined at ```framework/datastore/RelationsObject.h```:
      ```cpp
      template <class FROM> std::pair<FROM*, float> getRelatedFromWithWeight(const std::string& name = "",
      const std::string& namedRelation = "") const
      {
      RelationEntry entry = DataStore::Instance().getRelationWith(DataStore::c_FromSide, this, m_cacheDataStoreEntry, m_cacheArrayIndex, FROM::Class(), name, namedRelation);
      return std::make_pair(static_cast<FROM*>(entry.object), entry.weight);
      }
      ```
    hence ```getKlId()``` returns the second element of the pair which is just the weight. 


* ```m_KLMMCMom```- momentum of matched mc particle
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` : ```m_KLMMCMom = part->getMomentum().R();```
    ```getMomentum``` is defined at ```mdst/dataobjectsMCParticle.h``` as
    ```cpp
    ROOT::Math::XYZVector getMomentum() const
    {
      return ROOT::Math::XYZVector(m_momentum_x, m_momentum_y, m_momentum_z);
    }
    ```
<**Note** find out how they defined ```m_momentum_x, m_momentum_y, m_momentum_z``` as well as what ```.R()``` is. I suspect that it is a magnitude but not sure. 


* ```m_KLMMCPhi```- phi of matched mc particle
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMMCPhi        = part->getMomentum().Phi();```
    The ```getMomentum()``` is te same as above. ```Phi()``` is defined in ```framework/geometry/B2Vector3.h``` as ```DataType Phi() const { return X() == 0 && Y() == 0 ? 0 : atan2(Y(), X()); }```. This returns the azimuthal angle. 

* ```m_KLMMCTheta```- theta of matched mc particle
  - TYPE: ```Float_t```
  - This varibale is defined in ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMMCTheta      = part->getMomentum().Theta();``` where ```Theta()``` is defined in ```framework/geometry/B2Vector3.h``` as ``` DataType Theta() const { return X() == 0 && Y() == 0 && Z() == 0 ? 0 : atan2(Perp(), Z()); }```. This returns the MC polar angle. 
 
* ```m_KLMMom```- measured momentum
  - TYPE: ```Float_t```
>**Note** cannot find info about this 

* ```m_KLMPhi```- measured phi
  - TYPE: ```Float_t```
  - This variable is found at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMPhi = clusterPos.Phi();``` where ```clusterPos``` is defined as ```const ROOT::Math::XYZVector& clusterPos = cluster.getClusterPosition();``` where ```cluster``` is defined at the for loop:``` for (const KLMCluster& cluster : m_klmClusters) {``` 

* ```m_KLMTheta```- measured theta
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMTheta = clusterPos.Theta();``` where ```clusterPos``` is defined as ```const ROOT::Math::XYZVector& clusterPos = cluster.getClusterPosition();``` where ```cluster``` is defined at the for loop:``` for (const KLMCluster& cluster : m_klmClusters) {``` (same as above)
    ```Theta()``` is defined in ```framework/geometry/B2Vector3.h``` as ``` DataType Theta() const { return X() == 0 && Y() == 0 && Z() == 0 ? 0 : atan2(Perp(), Z()); }```. This returns the measured polar angle. 


* ```m_KLMMCStatus```- MC particles status
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMMCStatus     = part->getStatus();``` where ```getStatus``` is defined at ``` mdst/dataobjects/ECLCluster.h``` as: ```unsigned short getStatus() const {return m_status;}```. ```m_status``` is defined on the same code as ```void setStatus(EStatusBit status) { m_status = static_cast<short unsigned>(status);} ``` where the class ```EStatusBit``` is defined as:
      ```cpp
      /** bit 0: ECLCluster is matched to a ECL trigger cluster */
      c_TriggerCluster   = 1 << 0,
      /** bit 1: ECLCluster to ECLTRGCluster matcher was run */
      c_TriggerClusterMatching = 1 << 1,
      /** bit 2: ECLCluster has pulse shape discrimination variables.*/
      c_PulseShapeDiscrimination = 1 << 2,
      /** bit 3: ECLCluster has fit time that failed.*/
      c_fitTimeFailed = 1 << 3,
      /** bit 4: ECLCluster has time resolution calculation that failed.*/
      c_timeResolutionFailed = 1 << 4
      ```
    according to this, this is likely translated into an integer and fed into ```setStatus()```. 
    


* ```m_KLMMCLifetime```- MC partilces life time
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMMCLifetime   = part->getLifetime();``` where ```getLifetime()``` is defined at ```mdst/dataobjects/MCParticle.h``` as: ```float getLifetime() const { return m_decayTime - m_productionTime; }```. It notes that this returns the lifetime in ns, and if particle crosses the siulation boundary, it is set as the time spent in that volume.
    ```m_decayTime``` is defined as ```float getDecayTime() const { return m_decayTime; }``` where ``` void setDecayTime(float time) { m_decayTime = time; }```.
    ```m_productionTime``` is defined as ```float getProductionTime() const { return m_productionTime; }``` where ```  void setProductionTime(float time) { m_productionTime = time; }```
<**Note** NEXT STEP IS TO FIND OUT HOW THE VARIABLES ARE DEFINED

  


* ```m_KLMMCPDG```- pdg code of matched MCparticle
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMMCPDG        = std::abs(part->getPDG());``` where ```getPDG()``` is defined at ```mdst/dataobjects/MCParticle.h``` as ```int getPDG() const { return m_pdg; }```. ```m_pdg``` is defined as ```void setPDG(int pdg) { m_pdg = pdg; }```
 
<**Note** FIND OUT WHERE THE pdg INPUT IS

* ```m_KLMMCPrimaryPDG```- pdg code of MCparticles mother, for example pi0 for some gammas
  - TYPE: ```Float_t````
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMMCPrimaryPDG = getPrimaryPDG(part);```. ```getPrimaryPDG()``` is defined at ```reconstruction/modules/KlId/KLMExpert/KlId.h``` as a while loop which runs whilst True, checks whether Klong, Kshort, kaon, neutron, proton, pion, pi0, muon, electron or photon, then checks whether the particle has Mother Particle, if not, it stops the while loop and returns the pdg code of the particle. If there is a mother particle, it will do the above procedure again and return the pdg code of the mother particle. If the mother particle has a mother then it repeats etc. The code mentions that this is "very imprecise but suffcient to understand if the backgrounds are charged, hadronic or gammas which is whats relevant for KlId investigations.".
        


* ```m_KLMECLHypo```- hypotheis id of closest ecl cluster 5: gamma, 6:hadron
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMECLHypo = closestECLCluster->getHypotheses();```. ```getHypotheses()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ``` unsigned short getHypotheses() const {return m_hypotheses;}```. ```m_hypotheses``` is defined in the same code as the output of ```void setHypothesis(EHypothesisBit hypothesis) { m_hypotheses = static_cast<short unsigned>(hypothesis); }```. ```EHypothesisBit``` is a class defined in the same code as follows:
      ```cpp
      enum class EHypothesisBit : short unsigned int {
        /** None as initializer */
        c_none = 1 << 0,
        /** CR is split into a muon and n photons (T1) */
        c_muonNPhotons = 1 << 1,
        /** CR is reconstructed as a charged hadron (T2) */
        c_chargedHadron = 1 << 2,
        /** CR is split into an electron and n photons (T3) */
        c_electronNPhotons = 1 << 3,
        /** CR is split into n photons (N1) */
        c_nPhotons = 1 << 4,
        /** CR is reconstructed as a neutral hadron (N2) */
        c_neutralHadron = 1 << 5,
        /** CR is reconstructed as merged pi0 (N3)*/
        c_mergedPi0 = 1 << 6
      };
      ```
    




 
* ```m_KLMECLZMVA```- zernike mva output for closest ECL cluster (based on around 10 z-moments)
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMECLZMVA           = closestECLCluster->getZernikeMVA();```. ```getZernikeMVA()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ```double getZernikeMVA() const { return m_zernikeMVA; }``` which "Return MVA based hadron/photon value based on Zernike moments (shower shape variable)". ```m_zernikeMVA``` is defined as the output of ```void setZernikeMVA(double zernikemva) { m_zernikeMVA = zernikemva; }```
    In ```ecl/moduleseclFinalize/ECLFinalizerModule.cc```, it states ```eclCluster->setZernikeMVA(eclShower->getZernikeMVA());```. The corresponding ```setZernikeMVA()``` seem to be used in ```ecl/modules/eclShowerShape/ECLShowerShapeModule.cc``` to define m_zernikeMVA.

<**Note** has N1 N2 shower jargon I dont understand yet. Look into it when I have time
    
 
* ```m_KLMECLZ40```- zernike moment 4,0 of closest ecl cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMECLZ40            = closestECLCluster->getAbsZernike40();```. ```getZernike40()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ``` double getAbsZernike40() const { return m_absZernike40; }```, where ```m_absZernike40``` is defined as the output of ```setZernike40(double zernike40)```. This function is used at ```ecl/modules/eclFinalize/ECLFinalizerModule.cc```: ``` eclCluster->setAbsZernike40(eclShower->getAbsZernikeMoment(4, 0));```. Now looking for how ```getAbsZernikeMoment``` is defined in ```ecl/dataobjects/ECLShower.h```: ```  double getAbsZernikeMoment(unsigned int n, unsigned int m) const { return m_absZernikeMoments[(n * (n + 1)) / 2 + m - 1]; }``` where ```m_absZernikeMoments[(n * (n + 1)) / 2 + m - 1``` is set with the function ```setAbsZernikeMoment```: ```void setAbsZernikeMoment(unsigned int n, unsigned int m, double absZernikeMoment) { m_absZernikeMoments[(n * (n + 1)) / 2 + m - 1] = absZernikeMoment; }```. ```setAbsZernikeMoment``` is used in ```ecl/modules/eclShowerShape/ECLShowerShapeModule.cc``` as:
      ```cpp
      for (unsigned int n = 1; n <= 5; n++) {
        for (unsigned int m = 0; m <= n; m++) {
          eclShower->setAbsZernikeMoment(n, m, computeAbsZernikeMoment(projectedECLDigits, sumEnergies, n, m, rho0));
      }
      ```
    where ```computeAbsZernikeMoment``` is defined as:
      ```cpp
      std::complex<double> sum(0.0, 0.0);

      for (const auto projectedDigit : projectedDigits) {
        double normalizedRho = projectedDigit.rho / rho0;     // Normalize radial distance according to rho0.
        //Ignore crystals with rho > rho0, if requested
        if (normalizedRho > 1.0) {
          if (!m_zernike_useFarCrystals) continue;
          else normalizedRho = 1.0; //crystals with rho > rho0 are scaled to rho0 instead of discarded
        }

        sum += projectedDigit.energy * std::conj(zernikeValue(normalizedRho, projectedDigit.alpha, n, m));
      }
      return (n + 1.0) / TMath::Pi() * std::abs(sum) / totalEnergy;
      ```
>**Note** I have no idea what projectedDigit is.come back to it 

  


* ```m_KLMECLZ51```- zernike moment 5,1 of closest ECL cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as:```m_KLMECLZ51            = closestECLCluster->getAbsZernike51();``` where ```getAbsZernike51()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ```double getAbsZernike51() const { return m_absZernike51; }```. ```m_absZernike51``` is the output of ```void setAbsZernike51(double zernike51) { m_absZernike51 = zernike51; }```. Again, this function is used at ```ecl/modules/eclFinalize/ECLFinalizerModule.cc``` like the this ```eclCluster->setAbsZernike51(eclShower->getAbsZernikeMoment(5, 1));```. This goes back to ```getAbsZernikeMoment``` from above.
 
>**Note** COME BACK TO IT WHEN THE PREVIOUS FUNCTION IS WELL UNDERSTOOD AS IT USES THE SAME FUNCION THAT I CURRENTLY HAVE NO UNDERSTANDING OF

* ```m_KLMECLUncertaintyPhi```- phi uncertainty oof closeest ecl cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as: ```m_KLMECLUncertaintyPhi = closestECLCluster->getUncertaintyPhi();```. ```getUncertaintyPhi()``` is defined as follows in ```mdst/dataobjects/ECLCluster.h```: ```double getUncertaintyPhi() const { return (m_sqrtcovmat_11);}```. ```m_sqrtcovmat_11``` is defined in a function defined in the same code:
      ```cpp
      void setCovarianceMatrix(double covArray[6]){
      m_sqrtcovmat_11 = sqrt(fabs(covArray[2])); // phi
      }
      ```
    This function is used in ```ecl/modules/eclFinalize/ECLFinalizerModule.cc```:
      ```cpp
      double covmat[6] = {
        eclShower->getUncertaintyEnergy()* eclShower->getUncertaintyEnergy(),
        0.0,
        eclShower->getUncertaintyPhi()* eclShower->getUncertaintyPhi(),
        0.0,
        0.0,
        eclShower->getUncertaintyTheta()* eclShower->getUncertaintyTheta()
      };
      eclCluster->setCovarianceMatrix(covmat);
      ```
<**Note** continue later from looking for getUncertainty functions. perhaps I have already defined it earlier as it sounds familiar
      
 
* ```m_KLMECLUncertaintyTheta```- theta uncertainty of closest ECL cluster
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMECLUncertaintyTheta = closestECLCluster->getUncertaintyTheta();```. ```getUncertaintyTheta()``` is defined at ```mdst/dataobjects/ECLCluster.h``` as ``` double getUncertaintyTheta() const {return (m_sqrtcovmat_22);}```. the output is defined here:
      ```cpp
      void setCovarianceMatrix(double covArray[6]){
        m_sqrtcovmat_22 = sqrt(fabs(covArray[5])); // theta
      }
      ```
    where ```covArray``` is defined the same as above.
<**Note** WHAT IS ```fabs```
        


* ```m_KLMMCWeight```- mc weight
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as
       ```cpp
       const auto mcParticleWeightPair = cluster.getRelatedToWithWeight<MCParticle>();
       m_KLMMCWeight = mcParticleWeightPair.second;
       ```
      where ```getRelatedToWithWeight``` is defined at ```framework/datastore/RelationsObject.h``` as
      ```cpp
      template <class TO> std::pair<TO*, float> getRelatedToWithWeight(const std::string& name = "", const std::string& namedRelation = "") const{
        RelationEntry entry = DataStore::Instance().getRelationWith(DataStore::c_ToSide, this, m_cacheDataStoreEntry, m_cacheArrayIndex, TO::Class(), name, namedRelation);
        return std::make_pair(static_cast<TO*>(entry.object), entry.weight);
      }
      ```
      
<**Note** CONTINUE LATER
        


* ```m_KLMtrackFlag```- track flag for belle comparision
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMtrackFlag                   = cluster.getAssociatedTrackFlag();```, where ```getAssociatedTrackFlag()``` is defined at ```mdst/dataobjects/KLMCluster.h``` as
      ```cpp
      bool KLMCluster::getAssociatedTrackFlag() const{
        RelationVector<Track> tracks = this->getRelationsFrom<Track>();
        return tracks.size() > 0;
      }
      ```
    where ```getRelationsFrom``` is defined at ```framework/datastore/RelationsObject.h```:
      ```cpp
      template <class FROM> RelationVector<FROM> getRelationsFrom(const std::string& name = "", const std::string& namedRelation = "") const{
        return RelationVector<FROM>(DataStore::Instance().getRelationsWith(DataStore::c_FromSide, this, m_cacheDataStoreEntry, m_cacheArrayIndex, FROM::Class(), name, namedRelation));
      }
      ```
<**Note** CONTINUE ON LATER


      


* ```m_KLMeclFlag```- ecl flag for belle comparision
  - TYPE: ```Float_t```
  - This variable is defined at ```reconstruction/modules/KlId/DataWriter/DataWriterModule.cc``` as ```m_KLMeclFlag                     = cluster.getAssociatedEclClusterFlag();```. ```getAssociatedEclClusterFlag()``` is defined at ```mdst/dataobjects/KLMCluster.cc``` as
      ```cpp
      bool KLMCluster::getAssociatedEclClusterFlag() const{
        RelationVector<ECLCluster> eclClusters = this->getRelationsTo<ECLCluster>();
        return eclClusters.size() > 0;
      }
      ```
    where ```getRelationsTo``` is defined at ```framework/datastore/RelationsObject.h```:
      ```cpp
      template <class TO> RelationVector<TO> getRelationsTo(const std::string& name = "", const std::string& namedRelation = "") const{
        return RelationVector<TO>(DataStore::Instance().getRelationsWith(DataStore::c_ToSide, this, m_cacheDataStoreEntry, m_cacheArrayIndex, TO::Class(), name, namedRelation));
      }
      ```
>**Note** CONTINUE LATER


    

 

Other variables: 

* ```m_klmClusters```
  - TYPE: FIND IT
  - Defined at ```analysis/modules/RestOfEventBuilder/RestOfEventBuilderModule.cc```:
      ```cpp
      const KLMCluster* klmCluster = m\_klmClusters[i];
      ```
    where KLMCluster is a class defined at ```mdst/dataobjects/KLMCluster.h```


 



    



    
    
    

  

    
    
    









    




