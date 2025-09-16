# scripts/sbdb_client.py
from __future__ import annotations

import time
import requests
from typing import Dict, Any, Optional, List, Tuple, Sequence

_SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
_CAD  = "https://ssd-api.jpl.nasa.gov/cad.api"
_UA   = "neo-updater (+https://github.com/mcranny/neo-updater)"

DEFAULT_TIMEOUT = 25
DEFAULT_RETRIES = 2


# ------------------------- HTTP helper ------------------------- #

def _get_json(url: str, params: Dict[str, Any],
              timeout: int = DEFAULT_TIMEOUT,
              retries: int = DEFAULT_RETRIES) -> Dict[str, Any]:
    headers = {"User-Agent": _UA, "Accept": "application/json"}
    last: Optional[BaseException] = None
    for a in range(retries + 1):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception as e:  # noqa: BLE001 - surfacing last error
            last = e
            if a < retries:
                time.sleep(0.5 * (2 ** a))
            else:
                raise
    if last is not None:
        raise last
    raise RuntimeError("unreachable in _get_json")  # for static analyzers


# ------------------------- utils ------------------------- #

def _elements_list_to_map(elements: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for e in elements or []:
        name = e.get("name") or e.get("label")
        if name:
            out[str(name)] = e
    return out

def _f(x: Any) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


# ------------------------- SBDB (elements) ------------------------- #

def sbdb_lookup(*, sstr: str | None = None,
                des: str | None = None,
                spk: int | None = None,
                full_prec: bool = True,
                soln_epoch: bool = False,
                cd_epoch: bool = True,
                cd_tp: bool = True,
                phys_par: bool = True) -> Dict[str, Any]:
    """
    Raw SBDB call. One of sstr|des|spk is required.
    We default to full precision at the *MPC epoch* so we keep (a, ma).
    We also request cd-epoch and cd-tp to get human-readable times.
    """
    if not (sstr or des or spk):
        raise ValueError("sbdb_lookup: one of sstr|des|spk is required")
    p: Dict[str, Any] = {}
    if sstr: p["sstr"] = sstr
    if des:  p["des"]  = des
    if spk:  p["spk"]  = str(spk)
    if full_prec:  p["full-prec"]  = "true"
    if soln_epoch: p["soln-epoch"] = "true"
    if cd_epoch:   p["cd-epoch"]   = "true"
    if cd_tp:      p["cd-tp"]      = "true"
    if phys_par:   p["phys-par"]   = "true"
    return _get_json(_SBDB, p)


def normalize_sbdb_orbit(j: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize SBDB orbit into fields the planner/viewer need.
    Angles in deg, distances in AU, times in JD TDB.

    Output keys (subset, most important first):
      - a_AU, e, i_deg, raan_deg, argp_deg, ma_deg (may be None),
        tp_jd_tdb (if available), epoch_jd_tdb
      - cd_epoch (string), cd_tp (string if cd-tp requested)
      - period_d, period_y, n_deg_d
      - H (if present), orbit_id, equinox, source, producer
      - condition_code, moid_au, moid_jup_au, t_jup
    """
    orbit = j.get("orbit") or {}
    els_map = _elements_list_to_map(orbit.get("elements", []))

    def g(name: str) -> Optional[float]:
        d = els_map.get(name)
        return _f(d.get("value")) if d else None

    a    = g("a")
    e    = g("e")
    q    = g("q")
    Q    = g("ad")
    inc  = g("i")
    om   = g("om")
    w    = g("w")
    ma   = g("ma")
    per_d = g("per")
    n     = g("n")
    tp    = g("tp")  # TDB JD (float) when present

    epoch_jd_tdb = _f(orbit.get("epoch"))
    cd_epoch     = orbit.get("cd_epoch")

    # cd(t_p) appears as "cd_value" on the tp element when cd-tp=1
    cd_tp = None
    if "tp" in els_map:
        cd_tp = els_map["tp"].get("cd_value")

    # Physical H (if present)
    H_val = None
    phys = j.get("phys_par")
    if isinstance(phys, list):
        for p in phys:
            if p.get("name") == "H":
                H_val = _f(p.get("value"))
                break

    out = {
        "e": e, "a_AU": a, "q_AU": q, "Q_AU": Q,
        "i_deg": inc, "raan_deg": om, "argp_deg": w, "ma_deg": ma,
        "epoch_jd_tdb": epoch_jd_tdb, "cd_epoch": cd_epoch,
        "period_d": per_d, "period_y": (per_d/365.25 if per_d else None),
        "n_deg_d": n, "tp_jd_tdb": tp, "cd_tp": cd_tp,
        "H": H_val,
        "orbit_id": orbit.get("orbit_id"), "equinox": orbit.get("equinox"),
        "source": orbit.get("source"), "producer": orbit.get("producer"),
        "condition_code": orbit.get("condition_code"),
        "moid_au": _f(orbit.get("moid")), "moid_jup_au": _f(orbit.get("moid_jup")),
        "t_jup": _f(orbit.get("t_jup")),
    }
    return out


def fetch_elements_jpl1(identifier: str) -> Dict[str, Any]:
    """
    Fetch elements (full precision, MPC epoch) and normalize.
    If sstr is ambiguous (code=300), follow with 'des' using the first match.
    """
    j = sbdb_lookup(sstr=identifier, full_prec=True, soln_epoch=False,
                    cd_epoch=True, cd_tp=True, phys_par=True)
    if j.get("code") == 300 and j.get("list"):
        pdes = j["list"][0].get("pdes")
        if pdes:
            j = sbdb_lookup(des=pdes, full_prec=True, soln_epoch=False,
                            cd_epoch=True, cd_tp=True, phys_par=True)
    return normalize_sbdb_orbit(j)


# ---- Back-compat wrapper for legacy callers ----
def fetch_elements(sstr: str) -> Optional[Dict[str, Any]]:
    """
    Back-compatible with your *old* function signature.
    Returns keys that legacy code expects:
      a_AU, e, i_deg, raan_deg, argp_deg, M_deg, epoch_jd, source
    Note:
      - epoch_jd is a TDB JD
      - we still *also* compute tp and ma@epoch via fetch_elements_jpl1 for new code
    """
    el = fetch_elements_jpl1(sstr)
    # Require the classic 6 elements for the legacy return
    need = ("a_AU", "e", "i_deg", "raan_deg", "argp_deg", "ma_deg")
    if not el or any(el.get(k) is None for k in need):
        return None
    return {
        "a_AU": el["a_AU"],
        "e": el["e"],
        "i_deg": el["i_deg"],
        "raan_deg": el["raan_deg"],
        "argp_deg": el["argp_deg"],
        "M_deg": el["ma_deg"],                  # legacy alias
        "epoch_jd": el["epoch_jd_tdb"],         # legacy alias (TDB JD)
        "source": el.get("source", "JPL SBDB"),
        # Tip: if a caller *also* wants tp, they should call fetch_elements_jpl1()
    }


# ------------------------- CAD (close-approach) ------------------------- #

def cad_query(*,
              des: str | None = None,
              spk: int | None = None,
              body: str = "Earth",
              date_min: str | None = None,
              date_max: str | None = None,
              dist_max: str | None = None,
              dist_min: str | None = None,
              min_dist_min: str | None = None,
              min_dist_max: str | None = None,
              h_min: float | None = None,
              h_max: float | None = None,
              sort: str = "date",
              limit: int | None = None,
              limit_from: int | None = None,
              include_fullname: bool = False,
              include_diameter: bool = False) -> Dict[str, Any]:
    p: Dict[str, Any] = {"sort": sort}
    if des is not None:         p["des"] = des
    if spk is not None:         p["spk"] = str(spk)
    if body:                    p["body"] = body
    if date_min:                p["date-min"] = date_min
    if date_max:                p["date-max"] = date_max
    if dist_max:                p["dist-max"] = dist_max
    if dist_min:                p["dist-min"] = dist_min
    if min_dist_min:            p["min-dist-min"] = min_dist_min
    if min_dist_max:            p["min-dist-max"] = min_dist_max
    if h_min is not None:       p["h-min"] = str(h_min)
    if h_max is not None:       p["h-max"] = str(h_max)
    if limit is not None:       p["limit"] = str(limit)
    if limit_from is not None:  p["limit-from"] = str(limit_from)
    if include_fullname:        p["fullname"] = "true"
    if include_diameter:        p["diameter"] = "true"
    return _get_json(_CAD, p)


def normalize_cad_rows(j: Dict[str, Any]) -> List[Dict[str, Any]]:
    # 'fields' is a list of strings naming the columns
    fields: List[str] = [str(x) for x in j.get("fields", [])]
    idx: Dict[str, int] = {name: i for i, name in enumerate(fields)}

    # CAD 'data' is a list of rows; each row is a sequence of scalars/strings
    data: List[Sequence[Any]] = j.get("data", [])
    out: List[Dict[str, Any]] = []

    def get_val(row: Sequence[Any], name: str) -> Any:
        i = idx.get(name)
        if i is None or i < 0 or i >= len(row):
            return None
        return row[i]

    def to_float(x: Any) -> Optional[float]:
        if x is None:
            return None
        if isinstance(x, (int, float)):
            return float(x)
        if isinstance(x, str):
            try:
                return float(x.strip())
            except ValueError:
                return None
        return None

    for row in data:
        out.append({
            "des":            get_val(row, "des"),
            "orbit_id":       get_val(row, "orbit_id"),
            "jd_tdb":         to_float(get_val(row, "jd")),
            "cd_tdb":         get_val(row, "cd"),
            "dist_au":        to_float(get_val(row, "dist")),
            "dist_min_au":    to_float(get_val(row, "dist_min")),
            "dist_max_au":    to_float(get_val(row, "dist_max")),
            "v_rel_kms":      to_float(get_val(row, "v_rel")),
            "v_inf_kms":      to_float(get_val(row, "v_inf")),
            "t_sigma_fmt":    get_val(row, "t_sigma_f"),
            "h":              to_float(get_val(row, "h")) if "h" in idx else None,
            "diameter_km":    to_float(get_val(row, "diameter")) if "diameter" in idx else None,
            "diameter_sigma_km": to_float(get_val(row, "diameter_sigma")) if "diameter" in idx else None,
            "fullname":       get_val(row, "fullname") if "fullname" in idx else None,
            "body":           get_val(row, "body") if "body" in idx else "Earth",
        })

    return out


def fetch_close_approaches(identifier: str, *,
                           body: str = "Earth",
                           date_min: str | None = None,
                           date_max: str | None = None,
                           dist_max: str | None = None,
                           limit: int | None = 60,
                           include_fullname: bool = False,
                           include_diameter: bool = False) -> List[Dict[str, Any]]:
    j = cad_query(des=identifier, body=body, date_min=date_min, date_max=date_max,
                  dist_max=dist_max, limit=limit,
                  include_fullname=include_fullname, include_diameter=include_diameter)
    return normalize_cad_rows(j)


def discover_earth_close_approaches(*,
                                    date_min: str = "now",
                                    date_max: str = "+60",
                                    dist_max: str = "0.05",
                                    limit: int = 2000,
                                    attach_elements: bool = True) -> List[Dict[str, Any]]:
    """
    CNEOS-style feed: list Earth close approaches within window.
    For each unique 'des', keep the next upcoming CA. Optionally attach JPL 1 elements.
    """
    j = cad_query(date_min=date_min, date_max=date_max, dist_max=dist_max,
                  body="Earth", sort="date", limit=limit,
                  include_fullname=True, include_diameter=True)
    rows = normalize_cad_rows(j)
    by_des: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        des = r.get("des")
        if des and des not in by_des:
            by_des[des] = r
    out: List[Dict[str, Any]] = []
    for des, next_ca in by_des.items():
        el = None
        if attach_elements:
            try:
                el = fetch_elements_jpl1(des)
            except Exception:
                el = None
        out.append({"des": des, "next_ca": next_ca, "elements": el})
    return out
