# Megacryst melt evolution

Fractional crystallisation (FC) and digestion-fractional crystallisation (DFC) of olivine for megacryst parental melts.

The calculation is in `olivine-fc-dfc.ts`. It is the browser port of the SD2 spreadsheet (Greene et al.).

- `runFc` traces MgO–FeO–Ni and the equilibrium olivine along an FC path (Matzen DNi, T in °C).
- `runDfc` / `runDfcFromFc` grow an olivine rim while assimilating orthopyroxene (Pilbeam-style grain growth + AFC).
- `PARENTAL_DEFAULTS` is the reduced megacryst starting melt used in SD2 Fig. 11.
- `meltFromOlivineMargin` rebuilds the DFC starting melt from the FC olivine at a chosen F.

Used by the Hygrohodos workshop (Megacryst → Fractional Crystallization / Digestion-Fractional Crystallization).
