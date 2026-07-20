"""Precomputed low-dimensional universal odd-primary formulas."""

from __future__ import annotations

import base64
import json
import zlib

TensorTerm = tuple[tuple[int, ...], ...]
CatalogKey = tuple[int, int, int, bool]


_PRECOMPUTED_TERMS: dict[CatalogKey, dict[TensorTerm, int]] = {
    (3, 0, 1, True): {
        ((0, 2), (0, 1), (1, 2)): 2,
    },
    (3, 0, 2, True): {
        ((0, 1, 2), (0, 1, 2), (1, 2, 3)): 2,
        ((0, 2, 3), (0, 1, 3), (1, 2, 3)): 2,
        ((0, 1, 2), (0, 2, 3), (0, 1, 3)): 1,
        ((0, 1, 2), (1, 2, 3), (1, 2, 3)): 1,
    },
    (5, 0, 1, True): {
        ((0, 2), (0, 2), (0, 2), (0, 1), (1, 2)): 4,
        ((0, 2), (0, 1), (1, 2), (0, 1), (1, 2)): 1,
        ((0, 2), (0, 1), (1, 2), (0, 2), (0, 2)): 4,
        ((1, 2), (0, 1), (1, 2), (0, 2), (0, 1)): 1,
    },
}


_ENCODED_TERMS: dict[CatalogKey, str] = {
    (3, 1, 3, False): (
        "c-n=MOBR422t_ya1}RvJ;kf^8BvwZt%aOn9<6?=3Y)HtU3L!8G_BPYK9Y?K5Z8c|fXT`Mjl"
        "F`|U>EkD(C;v>jnJZD5Y#27TFK2$`lz!rw5^`+E!`?Htf^{AY1)ttUV7U9Fz_4TdSF+BdFJ"
        "EEioL@dqA0vlupDi$Ca^jygbB-1KJJ%DfxS-("
    ),
    (3, 1, 4, False): (
        "c-oDbTXyRp3`Hm8j~lRY=aQBAPvgc-D3Z<=<}Yz}WzgeD(C72{{I@-}vej*WZO?r>wzi$y&"
        "--(K{_}cmkN5XE`F^@cUf<Wq>)RUn`FuxypDvPX>%5+>t@C=iw$9Hgk419r?boxlw}0N<7R"
        "fc6<8vgBUptyh$D1ruvPV<N=GHc{uF2aY+2ejs_82L7{`t^n@+s?#WY490vd8_N?9t3fjk}U"
        "QF4sD<Xnx4<8DE=fb42jk(ah$<*tB(|Hj5!`&MP_ce5N6h%&Rz(dBsICucJurQ!m;4ydz`W)"
        "k`rho2hPZbN}@=DOeNoMDmk$o{Z77UMfKL?CpP*k@|M^Qp_K(OLNldk|*PiY>Wz$eMdG%WRr"
        "V`wHZB=tD=&bDk_=9y^~(SX5@V?wU9r#KVv=m9%spz=hxncDn5W53YO{oDO#85e2k#C`IwT{"
        "=A5*>V@}rT-fi47&&Mfc#U}M5BE8b?6Mt$hd7C{h`|Wg5LF>NL$!~so@<cMOvwn5SW=KW#wj"
        "}j9B8Mt|yM99II%71e=y0ZWRa7hY(g~9Q%{GJdlR!<$*tM(Vl1>y^k5k3Y_HOG@^C~B2Jx"
        "*&p>s*86Q;)UTlKk9*((Vu=o)xuUWc_`RF+cajCau%k)I85l6GU>o53+(En2b4}lXq(4vv-%"
        "7;gOKDnL+2+?7XR!o5fwvFJ5VhyxV(&UTh~8B}Z0dUeb!lrXV7D3L=sx7VQF=KUIgbiAC=#R"
        "^IOW=8T-rU@d5JuJbL_d!{zj0a=p^cI|O;iDX0afVZDWdIj_OQ{W6+ker=iEjd(iOfp33I("
        "^O;?Aq2#)@H-I3-XfIW`SvqbeJbqOZwY`l9KWKUP2wwjFjclk(`e-ml~cCis+Ih{TZR8WSs"
        "n?=jklWJE}ObjGMGb`d@keq_y!K3!Pth{<I|DJX@PBsm&^0(iqWcf@00{C!Lyks<dRR;+oi"
        "`9eJBaTAS5XamD<!>dmCB%_>>f(%P(vO|4BIBhoj{WJ{9rFiAPUr1~=>=lzKG%kw;`7T;OX"
        "{Up!xlog!wr0qPV(+xJy-j=-cQcEf><|p0W`<a1Ajlyn7(hc?!)gw>Tq*{Di$RwR2X-Z))("
        "UiioH_ZurTS))cTUn=`d^0jmfH^rkYah~Kv)Sel1=DW$`VXbjp=A"
    ),
}


def precomputed_terms(
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
) -> dict[TensorTerm, int] | None:
    """Return stored universal terms for a low-dimensional operation."""
    terms = _PRECOMPUTED_TERMS.get((p, r, source_degree, bockstein))
    if terms is None:
        encoded = _ENCODED_TERMS.get((p, r, source_degree, bockstein))
        if encoded is None:
            return None
        terms = _decode_terms(encoded)
    return dict(terms)


def _decode_terms(encoded: str) -> dict[TensorTerm, int]:
    payload = json.loads(zlib.decompress(base64.b85decode(encoded)).decode())
    return {
        tuple(tuple(factor) for factor in tensor): coefficient
        for tensor, coefficient in payload
    }
