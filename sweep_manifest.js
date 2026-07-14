// Costim antigen×arm off-tumor sweep — from COSTIM_OFFTOX_ANALYSIS_2026-07-14.csv
// Extensible: add antigen blocks + copy sweep/<TAA>_<arm>_<organ>.png as they render.
window.SWEEP = {
  default: { antigen: 'EGFR', arm: 'CD28', organ: 'large_int' },
  organs: ['tumor', 'large_int', 'small_int', 'skin', 'liver', 'spleen'],
  organLabel: { tumor: 'Tumour', large_int: 'Large intestine', small_int: 'Small intestine', skin: 'Skin', liver: 'Liver', spleen: 'Spleen' },
  armLabel: { NOCOSTIM: 'No costim', CD28: 'CD28', CD40: 'CD40', TNFRSF9: '4-1BB', CD27: 'CD27' },
  antigens: {
    EGFR: {
      label: 'EGFR', note: 'solid · baseline TI 1.29', baseline: 'NOCOSTIM',
      arms: ['NOCOSTIM', 'CD28', 'CD40', 'TNFRSF9', 'CD27'],
      arm: {
        NOCOSTIM: { off: 0, ti: 1.29 },
        CD28: { off: 57, ti: 0.89 },
        CD40: { off: 4, ti: 1.26 },
        TNFRSF9: { off: 8, ti: 1.22 },
        CD27: { off: 7, ti: 1.24 }
      }
    },
    CEACAM5: {
      label: 'CEACAM5', note: 'solid CRC · baseline TI 0.54 (below 1)', baseline: 'NOCOSTIM',
      arms: ['NOCOSTIM', 'CD28', 'CD40', 'TNFRSF9', 'CD27'],
      arm: {
        NOCOSTIM: { off: 0, ti: 0.54 },
        CD28: { off: 40, ti: 0.40 },
        CD40: { off: 3, ti: 0.53 },
        TNFRSF9: { off: 6, ti: 0.51 },
        CD27: { off: 4, ti: 0.52 }
      }
    },
    CEACAM6: {
      label: 'CEACAM6', note: 'solid CRC · baseline TI 0.59 (below 1)', baseline: 'NOCOSTIM',
      arms: ['NOCOSTIM', 'CD28', 'CD40', 'TNFRSF9', 'CD27'],
      arm: {
        NOCOSTIM: { off: 0, ti: 0.59 },
        CD28: { off: 48, ti: 0.41 },
        CD40: { off: 3, ti: 0.57 },
        TNFRSF9: { off: 7, ti: 0.55 },
        CD27: { off: 5, ti: 0.56 }
      }
    }
    // CD19 (heme reference) — incoming
  },
  incoming: ['CD19'],
  file: function (antigen, arm, organ) { return 'sweep/' + antigen + '_' + arm + '_' + organ + '.png'; }
};
