from __future__ import division
from enum import IntEnum


class ElementTags(IntEnum):
    TH_START = 0
    WM = 1
    GM = 2
    CSF = 3
    BONE = 4
    SCALP = 5
    EYE_BALLS = 6
    COMPACT_BONE = 7
    SPONGY_BONE = 8
    BLOOD = 9
    MUSCLE = 10
    EDEMA = 11
    NECROSIS = 16
    BURRHOLES = 19
    RESECTION_CAVITIES = 20
    TUMOR_SHELL = 24
    RESIDUAL_TUMOR = 25
    ELECTRODE_RUBBER_START = 100
    ELECTRODE_RUBBER = 100
    ELECTRODE_RUBBER_END = 499
    SALINE_START = 500
    SALINE = 500
    SALINE_END = 899
    TH_END = 999

    TH_SURFACE_START = 1000
    WM_TH_SURFACE = TH_SURFACE_START + WM
    GM_TH_SURFACE = TH_SURFACE_START + GM
    CSF_TH_SURFACE = TH_SURFACE_START + CSF
    BONE_TH_SURFACE = TH_SURFACE_START + BONE
    SCALP_TH_SURFACE = TH_SURFACE_START + SCALP
    EYE_BALLS_TH_SURFACE = TH_SURFACE_START + EYE_BALLS
    COMPACT_BONE_TH_SURFACE = TH_SURFACE_START + COMPACT_BONE
    SPONGY_BONE_TH_SURFACE = TH_SURFACE_START + SPONGY_BONE
    BLOOD_TH_SURFACE = TH_SURFACE_START + BLOOD
    MUSCLE_TH_SURFACE = TH_SURFACE_START + MUSCLE
    EDEMA_TH_SURFACE = TH_SURFACE_START + EDEMA
    NECROSIS_TH_SURFACE = TH_SURFACE_START + NECROSIS
    BURRHOLES_TH_SURFACE = TH_SURFACE_START + BURRHOLES
    RESECTION_CAVITIES_TH_SURFACE = TH_SURFACE_START + RESECTION_CAVITIES
    TUMOR_SHELL_TH_SURFACE = TH_SURFACE_START + TUMOR_SHELL
    RESIDUAL_TUMOR_TH_SURFACE = TH_SURFACE_START + RESIDUAL_TUMOR
    ELECTRODE_RUBBER_TH_SURFACE_START = TH_SURFACE_START + ELECTRODE_RUBBER_START
    ELECTRODE_RUBBER_TH_SURFACE = TH_SURFACE_START + ELECTRODE_RUBBER
    ELECTRODE_RUBBER_TH_SURFACE_END = TH_SURFACE_START + ELECTRODE_RUBBER_END
    SALINE_TH_SURFACE_START = TH_SURFACE_START + SALINE_START
    SALINE_TH_SURFACE = TH_SURFACE_START + SALINE
    SALINE_TH_SURFACE_END = TH_SURFACE_START + SALINE_END
    ELECTRODE_PLUG_SURFACE_START = 2000
    ELECTRODE_PLUG_SURFACE_END = 2499
    TH_SURFACE_END = 2499

    LH_SURFACE_START = 5000
    LH_WM_SURFACE = LH_SURFACE_START + WM
    LH_GM_SURFACE = LH_SURFACE_START + GM
    LH_SPHERE = LH_SURFACE_START + 100
    LH_SPHERE_REG = LH_SURFACE_START + 101
    LH_SURFACE_END = 5999

    LH_CENTRAL_SURFACE_START = 6000
    LH_CENTRAL_GM = LH_CENTRAL_SURFACE_START + GM
    LH_CENTRAL_LAYER_1 = LH_CENTRAL_SURFACE_START + 100
    LH_CENTRAL_LAYER_23 = LH_CENTRAL_SURFACE_START + 101
    LH_CENTRAL_LAYER_4 = LH_CENTRAL_SURFACE_START + 102
    LH_CENTRAL_LAYER_5 = LH_CENTRAL_SURFACE_START + 103
    LH_CENTRAL_LAYER_6 = LH_CENTRAL_SURFACE_START + 104
    LH_CENTRAL_SURFACE_END = 6999

    RH_SURFACE_START = 7000
    RH_WM_SURFACE = RH_SURFACE_START + WM
    RH_GM_SURFACE = RH_SURFACE_START + GM
    RH_SPHERE = RH_SURFACE_START + 100
    RH_SPHERE_REG = RH_SURFACE_START + 101
    RH_SURFACE_END = 7999

    RH_CENTRAL_SURFACE_START = 8000
    RH_CENTRAL_GM = RH_CENTRAL_SURFACE_START + GM
    RH_CENTRAL_LAYER_1 = RH_CENTRAL_SURFACE_START + 100
    RH_CENTRAL_LAYER_23 = RH_CENTRAL_SURFACE_START + 101
    RH_CENTRAL_LAYER_4 = RH_CENTRAL_SURFACE_START + 102
    RH_CENTRAL_LAYER_5 = RH_CENTRAL_SURFACE_START + 103
    RH_CENTRAL_LAYER_6 = RH_CENTRAL_SURFACE_START + 104
    RH_CENTRAL_SURFACE_END = 8999


class CentralLayerDepths:
    CENTRAL_LAYER_1 = 0.06
    CENTRAL_LAYER_23 = 0.4
    CENTRAL_LAYER_4 = 0.55
    CENTRAL_LAYER_5 = 0.65
    CENTRAL_LAYER_6 = 0.85


central_cortical_layer_depths: dict[int, float] = {
    ElementTags.LH_CENTRAL_LAYER_1: CentralLayerDepths.CENTRAL_LAYER_1,
    ElementTags.RH_CENTRAL_LAYER_1: CentralLayerDepths.CENTRAL_LAYER_1,
    ElementTags.LH_CENTRAL_LAYER_23: CentralLayerDepths.CENTRAL_LAYER_23,
    ElementTags.RH_CENTRAL_LAYER_23: CentralLayerDepths.CENTRAL_LAYER_23,
    ElementTags.LH_CENTRAL_LAYER_4: CentralLayerDepths.CENTRAL_LAYER_4,
    ElementTags.RH_CENTRAL_LAYER_4: CentralLayerDepths.CENTRAL_LAYER_4,
    ElementTags.LH_CENTRAL_LAYER_5: CentralLayerDepths.CENTRAL_LAYER_5,
    ElementTags.RH_CENTRAL_LAYER_5: CentralLayerDepths.CENTRAL_LAYER_5,
    ElementTags.LH_CENTRAL_LAYER_6: CentralLayerDepths.CENTRAL_LAYER_6,
    ElementTags.RH_CENTRAL_LAYER_6: CentralLayerDepths.CENTRAL_LAYER_6,
}

central_cortical_layer_names: dict[int, str] = {
    ElementTags.LH_CENTRAL_LAYER_1: "central_cl_1",
    ElementTags.RH_CENTRAL_LAYER_1: "central_cl_1",
    ElementTags.LH_CENTRAL_LAYER_23: "central_cl_23",
    ElementTags.RH_CENTRAL_LAYER_23: "central_cl_23",
    ElementTags.LH_CENTRAL_LAYER_4: "central_cl_4",
    ElementTags.RH_CENTRAL_LAYER_4: "central_cl_4",
    ElementTags.LH_CENTRAL_LAYER_5: "central_cl_5",
    ElementTags.RH_CENTRAL_LAYER_5: "central_cl_5",
    ElementTags.LH_CENTRAL_LAYER_6: "central_cl_6",
    ElementTags.RH_CENTRAL_LAYER_6: "central_cl_6",
}

central_cortical_layer_tags: list[int] = [
    ElementTags.LH_CENTRAL_LAYER_1,
    ElementTags.RH_CENTRAL_LAYER_1,
    ElementTags.LH_CENTRAL_LAYER_23,
    ElementTags.RH_CENTRAL_LAYER_23,
    ElementTags.LH_CENTRAL_LAYER_4,
    ElementTags.RH_CENTRAL_LAYER_4,
    ElementTags.LH_CENTRAL_LAYER_5,
    ElementTags.RH_CENTRAL_LAYER_5,
    ElementTags.LH_CENTRAL_LAYER_6,
    ElementTags.RH_CENTRAL_LAYER_6,
]

tissue_tags: list[int] = [
    ElementTags.WM,
    ElementTags.GM,
    ElementTags.CSF,
    ElementTags.BONE,
    ElementTags.SCALP,
    ElementTags.EYE_BALLS,
    ElementTags.COMPACT_BONE,
    ElementTags.SPONGY_BONE,
    ElementTags.BLOOD,
    ElementTags.MUSCLE,
    ElementTags.ELECTRODE_RUBBER,
    ElementTags.SALINE,
    ElementTags.EDEMA,
    ElementTags.NECROSIS,
    ElementTags.BURRHOLES,
    ElementTags.RESECTION_CAVITIES,
    ElementTags.TUMOR_SHELL,
    ElementTags.RESIDUAL_TUMOR
]

tissue_names: dict[int, str] = {
    ElementTags.WM: "WM",
    ElementTags.GM: "GM",
    ElementTags.CSF: "CSF",
    ElementTags.BONE: "Bone",
    ElementTags.SCALP: "Scalp",
    ElementTags.EYE_BALLS: "Eye_balls",
    ElementTags.COMPACT_BONE: "Compact_bone",
    ElementTags.SPONGY_BONE: "Spongy_bone",
    ElementTags.BLOOD: "Blood",
    ElementTags.MUSCLE: "Muscle",
    ElementTags.ELECTRODE_RUBBER: "Electrode_rubber",
    ElementTags.SALINE: "Saline",
    ElementTags.EDEMA: 'Edema',
    ElementTags.NECROSIS: 'Necrosis',
    ElementTags.BURRHOLES: 'Burrholes',
    ElementTags.RESECTION_CAVITIES: 'Resection Cavities',
    ElementTags.TUMOR_SHELL: 'Tumor shell',
    ElementTags.RESIDUAL_TUMOR: 'Residual Tumor'
}

tissue_conductivities: dict[int, float] = {
    ElementTags.WM: 0.126,
    ElementTags.GM: 0.275,
    ElementTags.CSF: 1.654,
    ElementTags.BONE: 0.010,
    ElementTags.SCALP: 0.465,
    ElementTags.EYE_BALLS: 0.5,
    ElementTags.COMPACT_BONE: 0.008,
    ElementTags.SPONGY_BONE: 0.025,
    ElementTags.BLOOD: 0.6,
    ElementTags.MUSCLE: 0.16,
    ElementTags.ELECTRODE_RUBBER: 29.4,
    ElementTags.SALINE: 1.0,
    ElementTags.EDEMA: 0.71,
    ElementTags.NECROSIS: 1.0,
    ElementTags.BURRHOLES: 1.654,
    ElementTags.RESECTION_CAVITIES: 1.654,
    ElementTags.TUMOR_SHELL: 0.24,
    ElementTags.RESIDUAL_TUMOR: 1.654
}

tissue_conductivity_descriptions: dict[int, str] = {
    ElementTags.WM: "brain white matter (from Wagner 2004)",
    ElementTags.GM: "brain gray matter (from Wagner 2004)",
    ElementTags.CSF: "cerebrospinal fluid (from Wagner 2004)",
    ElementTags.BONE: "average bone (from Wagner 2004)",
    ElementTags.SCALP: "average scalp (from Wagner 2004)",
    ElementTags.EYE_BALLS: "vitreous humour (from Opitz, Paulus, Thielscher, submitted)",
    ElementTags.COMPACT_BONE: "compact bone (from Opitz, Paulus, Thielscher, submitted)",
    ElementTags.SPONGY_BONE: "spongy bone (from Opitz, Paulus, Thielscher, submitted)",
    ElementTags.BLOOD: "Blood (from Gabriel et al, 2009)",
    ElementTags.MUSCLE: "Muscle (from Gabriel et al, 2009)",
    ElementTags.ELECTRODE_RUBBER: "for tDCS rubber electrodes",
    ElementTags.SALINE: "for tDCS sponge electrodes",
    ElementTags.EDEMA: 'Edema',
    ElementTags.NECROSIS: 'Necrosis',
    ElementTags.BURRHOLES: 'Burrholes (CSF)',
    ElementTags.RESECTION_CAVITIES: 'Resection Cavities (CSF)',
    ElementTags.TUMOR_SHELL: 'Tumor shell',
    ElementTags.RESIDUAL_TUMOR: 'Residual Tumor (CSF)'
}
