   switch(GetComputeIndex(isComputeProcess_dEdr,
                          isComputeProcess_d2Edr2,
                          isComputeEnergy,
                          isComputeForces,
                          isComputeParticleEnergy,
                          isComputeVirial,
                          isComputeParticleVirial))
   {
      case 0:
         ier = Compute< false, false,
                        false, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 1:
         ier = Compute< false, false,
                        false, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 2:
         ier = Compute< false, false,
                        false, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 3:
         ier = Compute< false, false,
                        false, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 4:
         ier = Compute< false, false,
                        false, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 5:
         ier = Compute< false, false,
                        false, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 6:
         ier = Compute< false, false,
                        false, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 7:
         ier = Compute< false, false,
                        false, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 8:
         ier = Compute< false, false,
                        false, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 9:
         ier = Compute< false, false,
                        false, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 10:
         ier = Compute< false, false,
                        false, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 11:
         ier = Compute< false, false,
                        false, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 12:
         ier = Compute< false, false,
                        false, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 13:
         ier = Compute< false, false,
                        false, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 14:
         ier = Compute< false, false,
                        false, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 15:
         ier = Compute< false, false,
                        false, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 16:
         ier = Compute< false, false,
                        true, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 17:
         ier = Compute< false, false,
                        true, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 18:
         ier = Compute< false, false,
                        true, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 19:
         ier = Compute< false, false,
                        true, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 20:
         ier = Compute< false, false,
                        true, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 21:
         ier = Compute< false, false,
                        true, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 22:
         ier = Compute< false, false,
                        true, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 23:
         ier = Compute< false, false,
                        true, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 24:
         ier = Compute< false, false,
                        true, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 25:
         ier = Compute< false, false,
                        true, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 26:
         ier = Compute< false, false,
                        true, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 27:
         ier = Compute< false, false,
                        true, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 28:
         ier = Compute< false, false,
                        true, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 29:
         ier = Compute< false, false,
                        true, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 30:
         ier = Compute< false, false,
                        true, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 31:
         ier = Compute< false, false,
                        true, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 32:
         ier = Compute< false, true,
                        false, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 33:
         ier = Compute< false, true,
                        false, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 34:
         ier = Compute< false, true,
                        false, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 35:
         ier = Compute< false, true,
                        false, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 36:
         ier = Compute< false, true,
                        false, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 37:
         ier = Compute< false, true,
                        false, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 38:
         ier = Compute< false, true,
                        false, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 39:
         ier = Compute< false, true,
                        false, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 40:
         ier = Compute< false, true,
                        false, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 41:
         ier = Compute< false, true,
                        false, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 42:
         ier = Compute< false, true,
                        false, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 43:
         ier = Compute< false, true,
                        false, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 44:
         ier = Compute< false, true,
                        false, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 45:
         ier = Compute< false, true,
                        false, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 46:
         ier = Compute< false, true,
                        false, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 47:
         ier = Compute< false, true,
                        false, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 48:
         ier = Compute< false, true,
                        true, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 49:
         ier = Compute< false, true,
                        true, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 50:
         ier = Compute< false, true,
                        true, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 51:
         ier = Compute< false, true,
                        true, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 52:
         ier = Compute< false, true,
                        true, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 53:
         ier = Compute< false, true,
                        true, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 54:
         ier = Compute< false, true,
                        true, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 55:
         ier = Compute< false, true,
                        true, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 56:
         ier = Compute< false, true,
                        true, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 57:
         ier = Compute< false, true,
                        true, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 58:
         ier = Compute< false, true,
                        true, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 59:
         ier = Compute< false, true,
                        true, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 60:
         ier = Compute< false, true,
                        true, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 61:
         ier = Compute< false, true,
                        true, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 62:
         ier = Compute< false, true,
                        true, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 63:
         ier = Compute< false, true,
                        true, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 64:
         ier = Compute< true, false,
                        false, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 65:
         ier = Compute< true, false,
                        false, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 66:
         ier = Compute< true, false,
                        false, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 67:
         ier = Compute< true, false,
                        false, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 68:
         ier = Compute< true, false,
                        false, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 69:
         ier = Compute< true, false,
                        false, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 70:
         ier = Compute< true, false,
                        false, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 71:
         ier = Compute< true, false,
                        false, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 72:
         ier = Compute< true, false,
                        false, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 73:
         ier = Compute< true, false,
                        false, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 74:
         ier = Compute< true, false,
                        false, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 75:
         ier = Compute< true, false,
                        false, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 76:
         ier = Compute< true, false,
                        false, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 77:
         ier = Compute< true, false,
                        false, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 78:
         ier = Compute< true, false,
                        false, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 79:
         ier = Compute< true, false,
                        false, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 80:
         ier = Compute< true, false,
                        true, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 81:
         ier = Compute< true, false,
                        true, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 82:
         ier = Compute< true, false,
                        true, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 83:
         ier = Compute< true, false,
                        true, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 84:
         ier = Compute< true, false,
                        true, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 85:
         ier = Compute< true, false,
                        true, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 86:
         ier = Compute< true, false,
                        true, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 87:
         ier = Compute< true, false,
                        true, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 88:
         ier = Compute< true, false,
                        true, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 89:
         ier = Compute< true, false,
                        true, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 90:
         ier = Compute< true, false,
                        true, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 91:
         ier = Compute< true, false,
                        true, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 92:
         ier = Compute< true, false,
                        true, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 93:
         ier = Compute< true, false,
                        true, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 94:
         ier = Compute< true, false,
                        true, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 95:
         ier = Compute< true, false,
                        true, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 96:
         ier = Compute< true, true,
                        false, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 97:
         ier = Compute< true, true,
                        false, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 98:
         ier = Compute< true, true,
                        false, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 99:
         ier = Compute< true, true,
                        false, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 100:
         ier = Compute< true, true,
                        false, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 101:
         ier = Compute< true, true,
                        false, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 102:
         ier = Compute< true, true,
                        false, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 103:
         ier = Compute< true, true,
                        false, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 104:
         ier = Compute< true, true,
                        false, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 105:
         ier = Compute< true, true,
                        false, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 106:
         ier = Compute< true, true,
                        false, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 107:
         ier = Compute< true, true,
                        false, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 108:
         ier = Compute< true, true,
                        false, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 109:
         ier = Compute< true, true,
                        false, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 110:
         ier = Compute< true, true,
                        false, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 111:
         ier = Compute< true, true,
                        false, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 112:
         ier = Compute< true, true,
                        true, false,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 113:
         ier = Compute< true, true,
                        true, false,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 114:
         ier = Compute< true, true,
                        true, false,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 115:
         ier = Compute< true, true,
                        true, false,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 116:
         ier = Compute< true, true,
                        true, false,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 117:
         ier = Compute< true, true,
                        true, false,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 118:
         ier = Compute< true, true,
                        true, false,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 119:
         ier = Compute< true, true,
                        true, false,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 120:
         ier = Compute< true, true,
                        true, true,
                        false, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 121:
         ier = Compute< true, true,
                        true, true,
                        false, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 122:
         ier = Compute< true, true,
                        true, true,
                        false, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 123:
         ier = Compute< true, true,
                        true, true,
                        false, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 124:
         ier = Compute< true, true,
                        true, true,
                        true, false,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 125:
         ier = Compute< true, true,
                        true, true,
                        true, false,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 126:
         ier = Compute< true, true,
                        true, true,
                        true, true,
                        false>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      case 127:
         ier = Compute< true, true,
                        true, true,
                        true, true,
                        true>(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy,
                  *virial,
                  particleVirial);
         break;
      default:
         std::cout << "Unknown compute function index" << std::endl;
         ier = true;
         break;
   }
