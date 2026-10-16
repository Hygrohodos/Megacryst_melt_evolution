/**
 * Browser port of artifacts/zz Internal/calculations/olivine_fc.py + olivine_dfc.py
 * (Greene et al. / SD2). Tracks MgO–FeO–Ni for the Mg#–Ni diagram.
 */

const MW_MGO = 40.3044;
const MW_FEO = 71.8444;
const MW_SIO2 = 60.08;

export type FcPoint = {
  f: number;
  t_c: number;
  mgo: number;
  feo: number;
  ni: number;
  ol_mg_num: number;
  ol_ni_ppm: number;
  dni: number;
  mgo_ol: number;
  feo_ol: number;
};

export type ParentalMelt = {
  kd: number;
  t0: number;
  dT: number;
  mgo0: number;
  feo0: number;
  ni0: number;
  maxF: number;
};

export type AssimilationMelt = {
  id: string;
  f: number;
  r: number;
  kd: number;
  dni: number;
};

/** Reduced megacryst FC used in SD2 Fig. 11. */
export const PARENTAL_DEFAULTS: ParentalMelt = {
  kd: 0.3,
  t0: 1264.93,
  dT: 1.649,
  mgo0: 28.5,
  feo0: 18.5065,
  ni0: 809.36,
  maxF: 30,
};

/** Reduced protokimberlite FC that feeds the DFC starting-melt sheet. */
export const PROTO_PARENTAL_DEFAULTS: ParentalMelt = {
  kd: 0.3,
  t0: 1276.91,
  dT: 8.0,
  mgo0: 26.5,
  feo0: 12.619,
  ni0: 690.67,
  maxF: 25,
};
export const ASSIM_DEFAULTS: Omit<AssimilationMelt, "id"> = {
  f: 7,
  r: 0.6,
  kd: 0.4,
  dni: 15,
};

export const MELT_COLORS = [
  "rgba(125, 155, 138, 0.95)",
  "#3b82f6",
  "#eab308",
  "#ef4444",
  "#a855f7",
  "#22c55e",
  "#f97316",
  "#06b6d4",
  "#ec4899",
  "#84cc16",
];

const ASSIM_OPX = { mgo: 34.0, feo: 5.9, ni: 700.0 };
export { ASSIM_OPX };

const MELT_ESTIMATE = {
  SiO2: 4.0,
  TiO2: 1.5,
  Al2O3: 0.4,
  FeO: 5.0,
  MgO: 28.0,
  CaO: 25.0,
  Na2O: 0.0,
  K2O: 2.0,
  H2O: 4.0,
  CO2: 30.0,
};
const RMM = {
  SiO2: 60.08,
  TiO2: 79.866,
  Al2O3: 101.961,
  FeO: 71.8444,
  MgO: 40.3044,
  CaO: 56.0774,
  Na2O: 61.9789,
  K2O: 94.196,
  H2O: 18.01528,
  CO2: 44.0095,
};

/** SD2 DFC starting melt: Fe–Mg equilibrium with FC olivine at F, not the FC bulk liquid. */
export function meltFromOlivineMargin(
  olMgo: number,
  olFeo: number,
  olNi: number,
  kd: number,
  dni: number,
): { mgo: number; feo: number; ni: number } {
  const olMgFe = olMgo / RMM.MgO / (olFeo / RMM.FeO);
  const meltMgFe = kd * olMgFe;
  const meltMgNum = meltMgFe / (1 + meltMgFe);
  const ar = MELT_ESTIMATE.FeO / RMM.FeO;
  const as_ = MELT_ESTIMATE.MgO / RMM.MgO;
  const ba = ar + as_;
  const feMol = ba * (1 - meltMgNum);
  const mgMol = ba * meltMgNum;
  const oxides = {
    SiO2: (MELT_ESTIMATE.SiO2 / RMM.SiO2) * RMM.SiO2,
    TiO2: (MELT_ESTIMATE.TiO2 / RMM.TiO2) * RMM.TiO2,
    Al2O3: ((MELT_ESTIMATE.Al2O3 * 2) / RMM.Al2O3) * 0.5 * RMM.Al2O3,
    FeO: feMol * RMM.FeO,
    MgO: mgMol * RMM.MgO,
    CaO: (MELT_ESTIMATE.CaO / RMM.CaO) * RMM.CaO,
    Na2O: ((2 * MELT_ESTIMATE.Na2O) / RMM.Na2O) * 0.5 * RMM.Na2O,
    K2O: ((2 * MELT_ESTIMATE.K2O) / RMM.K2O) * 0.5 * RMM.K2O,
    H2O: ((2 * MELT_ESTIMATE.H2O) / RMM.H2O) * 0.5 * RMM.H2O,
    CO2: (MELT_ESTIMATE.CO2 / RMM.CO2) * RMM.CO2,
  };
  const total = Object.values(oxides).reduce((a, b) => a + b, 0);
  return {
    mgo: (oxides.MgO * 100) / total,
    feo: (oxides.FeO * 100) / total,
    ni: olNi / dni,
  };
}

function olivineFromLiquid(kd: number, mgoL: number, feoL: number) {
  const mgFeL = mgoL / MW_MGO / (feoL / MW_FEO);
  const mgFeS = mgFeL / kd;
  const mgNum = mgFeS / (mgFeS + 1);
  const sio2 = MW_SIO2;
  const mgo = 2 * mgNum * MW_MGO;
  const feo = 2 * (1 - mgNum) * MW_FEO;
  const sum = sio2 + mgo + feo;
  return { mg_number_pct: mgNum * 100, mgo_s: (mgo / sum) * 100, feo_s: (feo / sum) * 100 };
}

function matzenDni(tC: number, mgoL: number, mgoS: number): number {
  if (tC <= 0 || mgoL <= 0 || mgoS <= 0) return NaN;
  return Math.exp(4338 / tC - 1.956 - Math.log(mgoL / (mgoS / 2)));
}

/** Olivine Ni (ppm) in Fe–Mg–Ni equilibrium with a given melt. Same Matzen DNi as the FC path. */
export function olivineNiPpm(mgo: number, feo: number, niPpm: number, tC: number, kd = 0.3): number {
  const ol = olivineFromLiquid(kd, mgo, feo);
  const dni = matzenDni(tC, mgo, ol.mgo_s);
  if (!Number.isFinite(dni) || !Number.isFinite(niPpm)) return NaN;
  return niPpm * dni;
}

export type TempRedox = "oxidized" | "reduced1" | "reduced2";
export type NiRedox = "oxidized" | "reducedA" | "reducedB" | "reduced";

/** SD2 / Greene et al. (2023) eqs. 3–5. Reduced option 1 is the oxidised fit minus 450 °C. */
export function temperatureFromMgo(mgo: number, redox: TempRedox): number {
  const oxidized = (mgo + 10.4683) / 0.022723;
  if (redox === "oxidized") return oxidized;
  if (redox === "reduced1") return oxidized - 450;
  return 12.9 * mgo + 899;
}

/** SD2 / Greene et al. (2023) eqs. 7–9. `reduced` is the mean of the two reduced fits. */
export function niFromFeo(feo: number, redox: NiRedox): number {
  if (redox === "oxidized") return 34.6 * feo + 83.6;
  const a = 48.9 * feo - 221;
  const b = 16.9 * feo + 622;
  if (redox === "reducedA") return a;
  if (redox === "reducedB") return b;
  return 0.5 * (a + b);
}

/** Starting melt at MgO/FeO = 1.54 and the olivine in Fe–Mg–Ni equilibrium with it. */
export function startingOlivine(
  mgo: number,
  redox: TempRedox,
  kd = 0.3,
): { mgo: number; feo: number; ni: number; t: number; olMg: number; olNi: number } {
  const feo = mgo / 1.54;
  const ni = niFromFeo(feo, redox === "oxidized" ? "oxidized" : "reduced");
  const t = temperatureFromMgo(mgo, redox);
  const ol = olivineFromLiquid(kd, mgo, feo);
  const dni = matzenDni(t, mgo, ol.mgo_s);
  return { mgo, feo, ni, t, olMg: ol.mg_number_pct, olNi: ni * dni };
}

export function runFc(p: ParentalMelt): FcPoint[] {
  const step = 1;
  const n = Math.floor(Math.max(0, p.maxF) / step) + 1;
  const out: FcPoint[] = [];
  let mgo = p.mgo0;
  let feo = p.feo0;
  let ni = p.ni0;
  let t = p.t0;
  for (let i = 0; i < n; i++) {
    const f = i * step;
    const ol = olivineFromLiquid(p.kd, mgo, feo);
    const dni = matzenDni(t, mgo, ol.mgo_s);
    const niOl = ni * dni;
    out.push({
      f,
      t_c: t,
      mgo,
      feo,
      ni,
      ol_mg_num: ol.mg_number_pct,
      ol_ni_ppm: niOl,
      dni,
      mgo_ol: ol.mgo_s,
      feo_ol: ol.feo_s,
    });
    if (i === n - 1) break;
    const dF = step / 100;
    mgo -= dF * ol.mgo_s;
    feo -= dF * ol.feo_s;
    ni -= dF * niOl;
    t -= p.dT * step;
  }
  return out;
}

export type DfcMelt = {
  mgo: number;
  feo: number;
  ni: number;
  kd: number;
  dni: number;
  r: number;
  assimMgo?: number;
  assimFeo?: number;
  assimNi?: number;
};

/** DFC path from an explicit starting melt and assimilant (SD2 DFC input). */
export function runDfc(melt: DfcMelt): { ol_mg_num: number; ol_ni_ppm: number }[] {
  const assim = {
    mgo: melt.assimMgo ?? ASSIM_OPX.mgo,
    feo: melt.assimFeo ?? ASSIM_OPX.feo,
    ni: melt.assimNi ?? ASSIM_OPX.ni,
  };
  const radiusCoreUm = 1200;
  const dxUm = 0.3;
  const maxRimUm = 80;
  const density = 3220;
  const rOx = 0.08;
  const rCO2 = 0.297;
  const fSp = 0.5;
  const wtpctOl = 34;
  const wtpctOpx = 6;
  const radiusCoreM = radiusCoreUm * 1e-6;
  const vOlC = (4 / 3) * Math.PI * radiusCoreM ** 3;
  const mOlC = vOlC * density;
  const mBi = mOlC * (100 / wtpctOl);
  let mMelt = mBi * ((100 - wtpctOl - wtpctOpx) / 100);
  let mOpx = mBi * (wtpctOpx / 100);
  const mMeltI = mMelt;
  let radiusM = radiusCoreM;
  let mOl = mOlC;
  let mgoL = melt.mgo;
  let feoL = melt.feo;
  let niL = melt.ni;
  const dxM = dxUm * 1e-6;
  const pts: { ol_mg_num: number; ol_ni_ppm: number }[] = [];
  const ol0 = olivineFromLiquid(melt.kd, mgoL, feoL);
  pts.push({ ol_mg_num: ol0.mg_number_pct, ol_ni_ppm: niL * melt.dni });
  let fOld = 1;
  const nSteps = Math.ceil(maxRimUm / dxUm) + 1;
  for (let step = 1; step <= nSteps; step++) {
    radiusM += dxM;
    if (radiusM - radiusCoreM > maxRimUm * 1e-6 + 1e-12) break;
    const vOl = (4 / 3) * Math.PI * radiusM ** 3;
    const mOlNew = vOl * density;
    const dmOl = mOlNew - mOl;
    mOl = mOlNew;
    const dmOx = rOx * dmOl;
    const dmLiq = (1 - fSp) * dmOx + fSp * dmOx + rCO2 * dmOl + dmOl;
    const dmAssim = dmLiq * melt.r;
    mMelt = mMelt + dmAssim - dmLiq;
    mOpx -= dmAssim;
    if (mMelt <= 0 || mOpx < 0) break;
    const fNew = mMelt / mMeltI;
    const olPrev = olivineFromLiquid(melt.kd, mgoL, feoL);
    const niS = niL * melt.dni;
    const afc = (cOld: number, cA: number, cS: number) => {
      if (Math.abs(melt.r - 1) < 1e-12 || fOld <= 0 || fNew <= 0) return cOld;
      return cOld + ((melt.r * cA - melt.r * cOld + cOld - cS) / (melt.r - 1)) * Math.log(fNew / fOld);
    };
    mgoL = afc(mgoL, assim.mgo, olPrev.mgo_s);
    feoL = afc(feoL, assim.feo, olPrev.feo_s);
    niL = afc(niL, assim.ni, niS);
    fOld = fNew;
    const ol = olivineFromLiquid(melt.kd, mgoL, feoL);
    pts.push({ ol_mg_num: ol.mg_number_pct, ol_ni_ppm: niL * melt.dni });
  }
  return pts;
}

function fcAtF(fc: FcPoint[], f: number): FcPoint | null {
  if (!fc.length) return null;
  const exact = fc.find((p) => Math.abs(p.f - f) < 1e-9);
  if (exact) return exact;
  let lo = fc[0]!;
  let hi = fc[fc.length - 1]!;
  for (const p of fc) {
    if (p.f <= f && p.f >= lo.f) lo = p;
    if (p.f >= f && p.f <= hi.f) hi = p;
  }
  if (Math.abs(hi.f - lo.f) < 1e-9) return lo;
  const t = (f - lo.f) / (hi.f - lo.f);
  const mix = (a: number, b: number) => a + (b - a) * t;
  return {
    f,
    t_c: mix(lo.t_c, hi.t_c),
    mgo: mix(lo.mgo, hi.mgo),
    feo: mix(lo.feo, hi.feo),
    ni: mix(lo.ni, hi.ni),
    ol_mg_num: mix(lo.ol_mg_num, hi.ol_mg_num),
    ol_ni_ppm: mix(lo.ol_ni_ppm, hi.ol_ni_ppm),
    dni: mix(lo.dni, hi.dni),
    mgo_ol: mix(lo.mgo_ol, hi.mgo_ol),
    feo_ol: mix(lo.feo_ol, hi.feo_ol),
  };
}

export function runDfcFromFc(
  fc: FcPoint[],
  melt: AssimilationMelt,
): { ol_mg_num: number; ol_ni_ppm: number }[] {
  const start = fcAtF(fc, melt.f);
  if (!start) return [];
  const recon = meltFromOlivineMargin(start.mgo_ol, start.feo_ol, start.ol_ni_ppm, melt.kd, melt.dni);
  return runDfc({
    mgo: recon.mgo,
    feo: recon.feo,
    ni: recon.ni,
    kd: melt.kd,
    dni: melt.dni,
    r: melt.r,
  });
}

export function parentalSummary(p: ParentalMelt) {
  const fc = runFc({ ...p, maxF: 0 });
  const row = fc[0];
  return {
    ol_mg_num: row?.ol_mg_num ?? NaN,
    ol_ni_ppm: row?.ol_ni_ppm ?? NaN,
    dni: row?.dni ?? NaN,
  };
}

const MW_NIO = 74.6928;
const MW_NI = 58.6934;

/** Bulk parental melt oxides: panel MgO–FeO–Ni plus the SD2 spreadsheet remainder. */
export function parentalBulkOxides(p: ParentalMelt): Record<string, number> {
  const mg = p.mgo0 / MW_MGO;
  const fe = p.feo0 / MW_FEO;
  const den = mg + fe;
  return {
    SiO2: MELT_ESTIMATE.SiO2,
    TiO2: MELT_ESTIMATE.TiO2,
    Al2O3: MELT_ESTIMATE.Al2O3,
    Cr2O3: 0,
    FeO: p.feo0,
    MnO: 0,
    MgO: p.mgo0,
    CaO: MELT_ESTIMATE.CaO,
    Na2O: MELT_ESTIMATE.Na2O,
    NiO: (p.ni0 * MW_NIO) / MW_NI / 1e4,
    "Mg#": den > 0 ? (100 * mg) / den : NaN,
  };
}

export function newMeltId(): string {
  return `m-${Math.random().toString(36).slice(2, 8)}`;
}

export function defaultMelt(_index: number, prev?: AssimilationMelt): AssimilationMelt {
  return {
    id: newMeltId(),
    f: prev ? Math.min(30, prev.f + 3) : ASSIM_DEFAULTS.f,
    r: prev?.r ?? ASSIM_DEFAULTS.r,
    kd: prev?.kd ?? ASSIM_DEFAULTS.kd,
    dni: prev?.dni ?? ASSIM_DEFAULTS.dni,
  };
}
