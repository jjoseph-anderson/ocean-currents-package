import math
from datetime import datetime
from typing import List, Tuple

def phase_velocity(period: float, direction: float) -> Tuple[float, float, float]:
    """Calculate phase velocity and its components from wave period and direction."""
    g = 9.81
    cph = period * g / (2 * math.pi)
    vx = -cph * math.sin(math.radians(direction))
    vy = -cph * math.cos(math.radians(direction))
    return cph, vx, vy

def current_speed_projections(speed: float, direction: float) -> Tuple[float, float]:
    """Project surface current speed into x and y components."""
    vx = -speed * math.sin(math.radians(direction))
    vy = -speed * math.cos(math.radians(direction))
    return vx, vy

def steepness_hl(height: float, length: int) -> float:
    """Calculate wave steepness as height over wavelength."""
    return height / length if length > 0 else 0.0

class Para:
    """Wave and current parameters from WAMOS."""
    lat: float
    lon: float
    dt: datetime
    year: int
    month: int
    day: int
    hour: int
    minutes: int
    seconds: int

    h: float
    tp: float
    tm: float
    lp: int
    dm: int
    dp: int
    vp: float
    vxp: float
    vyp: float
    vm: float
    vxm: float
    vym: float
    stp: float

    ps: float
    ds: int
    ls: int
    vs: float
    vxs: float
    vys: float

    pw: float
    dw: int
    lw: int
    vv: float
    vx: float
    vy: float

    usp: float
    dir: float
    vxu: float
    vyu: float

    hmax: float
    tlim: float

class Peak:
    """Peak sea and swell system parameters from WAMOS."""
    lat: float
    lon: float
    dt: datetime

    h: float
    tp: float
    dp: int
    lp: int
    vp: float
    vxp: float
    vyp: float
    stp: float

    hw: float
    pw: float
    dw: int
    lw: int
    vv: float
    vx: float
    vy: float
    sw: float

    hs1: float
    ps1: float
    ds1: int
    ls1: int
    vs1: float
    vxs1: float
    vys1: float
    ss1: float

    hs2: float
    ps2: float
    ds2: int
    ls2: int
    vs2: float
    vxs2: float
    vys2: float
    ss2: float

    hs3: float
    ps3: float
    ds3: int
    ls3: int
    vs3: float
    vxs3: float
    vys3: float
    ss3: float

    usp: float
    dir: float
    vxu: float
    vyu: float

    hmax: float
    tlim: float

# ─────────────────────────────────────────────────────────────
# Parsing Functions
# ─────────────────────────────────────────────────────────────

def split_record(line: str) -> Para:
    """Parse a line of WAMOS data into a Para object."""
    lat0, lon0 = 32.1306, 34.7872
    pack = line.split()
    ret = Para()

    dt_str = pack[0]
    ret.year, ret.month, ret.day = int(dt_str[:4]), int(dt_str[4:6]), int(dt_str[6:8])
    ret.hour, ret.minutes, ret.seconds = int(dt_str[8:10]), int(dt_str[10:12]), int(dt_str[12:14])
    ret.dt = datetime(ret.year, ret.month, ret.day, ret.hour, ret.minutes, ret.seconds)
    ret.lat, ret.lon = lat0, lon0

    ret.h = float(pack[1]) * 0.33 - 0.1
    ret.tp, ret.tm = float(pack[2]), float(pack[3])
    ret.lp, ret.dm, ret.dp = int(pack[4]), int(pack[5]), int(pack[6])
    ret.vp, ret.vxp, ret.vyp = phase_velocity(ret.tp, ret.dp)
    ret.vm, ret.vxm, ret.vym = phase_velocity(ret.tm, ret.dm)
    ret.stp = steepness_hl(ret.h, ret.lp)

    ret.ps, ret.ds, ret.ls = float(pack[7]), int(pack[8]), int(pack[9])
    ret.vs, ret.vxs, ret.vys = phase_velocity(ret.ps, ret.ds)

    ret.pw, ret.dw, ret.lw = float(pack[10]), int(pack[11]), int(pack[12])
    ret.vv, ret.vx, ret.vy = phase_velocity(ret.pw, ret.dw)

    ret.usp, ret.dir = float(pack[13]), float(pack[14])
    ret.vxu, ret.vyu = current_speed_projections(ret.usp, ret.dir)

    ret.hmax = float(pack[18]) * 0.33 - 0.1
    ret.tlim = float(pack[19])

    return ret

def read_para(path: str) -> List[Para]:
    """Read a WAMOS file and return a list of Para objects."""
    with open(path, 'r', encoding="utf8", errors="ignore") as fi:
        fi.readline()  # skip header
        return [split_record(line) for line in fi if line.strip()]

def split_peak(line: str) -> Peak:
    """Parse a line of peak WAMOS data into a Peak object."""
    lat0, lon0 = 32.1306, 34.7872
    pack = line.split()
    ret = Peak()

    dt_str = pack[0]
    year, month, day = int(dt_str[:4]), int(dt_str[4:6]), int(dt_str[6:8])
    hour, minutes, seconds = int(dt_str[8:10]), int(dt_str[10:12]), int(dt_str[12:14])
    ret.dt = datetime(year, month, day, hour, minutes, seconds)
    ret.lat, ret.lon = lat0, lon0

    ret.h = float(pack[1]) * 0.33 - 0.1
    ret.tp, ret.dp, ret.lp = float(pack[2]), int(pack[3]), int(pack[4])
    ret.vp, ret.vxp, ret.vyp = phase_velocity(ret.tp, ret.dp)
    ret.stp = steepness_hl(ret.h, ret.lp)

    ret.hw = float(pack[5]) * 0.33 - 0.1
    ret.pw, ret.dw, ret.lw = float(pack[6]), int(pack[7]), int(pack[8])
    ret.vv, ret.vx, ret.vy = phase_velocity(ret.pw, ret.dw)
    ret.sw = steepness_hl(ret.hw, ret.lw)

    for i in range(3):
        hs, ps, ds, ls = float(pack[9 + i*4]) * 0.33 - 0.1, float(pack[10 + i*4]), int(pack[11 + i*4]), int(pack[12 + i*4])
        vs, vxs, vys = phase_velocity(ps, ds)
        ss = steepness_hl(hs, ls)
        setattr(ret, f'hs{i+1}', hs)
        setattr(ret, f'ps{i+1}', ps)
        setattr(ret, f'ds{i+1}', ds)
        setattr(ret, f'ls{i+1}', ls)
        setattr(ret, f'vs{i+1}', vs)
        setattr(ret, f'vxs{i+1}', vxs)
        setattr(ret, f'vys{i+1}', vys)
        setattr(ret, f'ss{i+1}', ss)

    ret.usp, ret.dir = float(pack[21]), float(pack[22])
    ret.vxu, ret.vyu = current_speed_projections(ret.usp, ret.dir)
    ret.tlim = float(pack[25])

    return ret

def read_peak(path: str) -> List[Peak]:
    """Read a WAMOS peak file and return a list of Peak objects."""
    with open(path, 'r', encoding="utf8", errors="ignore") as fi:
        fi.readline()  # skip header
        return [split_peak(line) for line in fi if line.strip()]
