#ifndef __VAANGO_DEFORMATION_GRADIENT_COMPUTER_H__
#define __VAANGO_DEFORMATION_GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>

#include <vector>

namespace Uintah {
    class Patch;
    class DataWarehouse;
}

namespace Vaango {

    class PeridynamicsFlags;
    class PeridynamicsMaterial;
    class PeridynamicsLabel;

    class DeformationGradientComputer {
        public:

            DeformationGradientComputer(PeridynamicsFlags* flags, PeridynamicsLabel* labels);
            virtual ~DeformationGradientComputer();

            void addInitialComputesAndRequires(Uintah::Task* task,
                                               const PeridynamicsMaterial* matl,
                                               const Uintah::PatchSet*);

            void addComputesAndRequires(Uintah::Task* task,
                                        const PeridynamicsMaterial* matl,
                                        const Uintah::PatchSet*);


            void initializeGradient(PeridynamicsMaterial* matl, 
                                    const Uintah::Patch* patch,    
                                    Uintah::DataWarehouse* new_dw);

            void computeDeformationGradient(const Uintah::Patch* patch,
                                            const PeridynamicsMaterial* matl,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw);

        private:

            PeridynamicsLabel* d_labels;
            PeridynamicsFlags* d_flags;

            void addComputesAndRequiresExplicit(Uintah::Task* task,
                                                const PeridynamicsMaterial* matl);

            void initializeGradientImplicit(const Uintah::Patch* patch,
                                            const PeridynamicsMaterial* matl,
                                            Uintah::DataWarehouse* new_dw);

            void initializeGradientExplicit(const Uintah::Patch* patch,
                                            const PeridynamicsMaterial* matl,
                                            Uintah::DataWarehouse* new_dw);

            void computeDeformationGradientImplicit(const Uintah::Patch* patch,
                                                    const PeridynamicsMaterial* matl,
                                                    Uintah::DataWarehouse* old_dw,
                                                    Uintah::DataWarehouse* new_dw);

            void computeDeformationGradientExplicit(const Uintah::Patch* patch,
                                                    const PeridynamicsMaterial* matl,
                                                    Uintah::DataWarehouse* old_dw,
                                                    Uintah::DataWarehouse* new_dw);

            static const Uintah::Matrix3 Identity;
    };
}

#endif
