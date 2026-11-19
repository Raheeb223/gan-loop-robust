"""Build a 4-layer .kicad_pcb from a primitives JSON (scripts/kicad_prims.py).

Run with KiCad's bundled interpreter:
    "C:/Program Files/KiCad/10.0/bin/python.exe" scripts/kicad_build.py prim.json out.kicad_pcb
Solver coordinates (y up) are mapped to KiCad (y down) with an origin offset.
The physical stackup (copper and dielectric thicknesses) is written into the
board file's (setup (stackup ...)) section so that Board Setup shows it and
fabrication notes are generated from it.
"""
import json
import re
import sys

import pcbnew as P

OX, OY = 100.0, 100.0
LAYER = {"F.Cu": P.F_Cu, "In1.Cu": P.In1_Cu, "In2.Cu": P.In2_Cu, "B.Cu": P.B_Cu}


def V(x, y):
    return P.VECTOR2I(P.FromMM(OX + x), P.FromMM(OY - y))


def main(src, dst):
    d = json.load(open(src))
    b = P.BOARD()
    b.SetCopperLayerCount(4)
    ds = b.GetDesignSettings()
    st = d["stackup"]
    ds.SetBoardThickness(P.FromMM(st["T"] + st["cu"]))   # outer copper centre-to-centre + copper
    nets = {}
    names = ["VIN", "SW", "GND"] + sorted({pd["net"] for f in d["fets"] for pd in f["pads"]} - {"VIN", "SW", "GND"})
    for n in names:
        ni = P.NETINFO_ITEM(b, n)
        b.Add(ni)
        nets[n] = ni
    zones = []
    for r in d["rects"]:
        w, h = r["x1"] - r["x0"], r["y1"] - r["y0"]
        if min(w, h) < 0.3:
            # narrow copper (EPC finger) -> track along its long axis, width = short side
            t = P.PCB_TRACK(b)
            if w >= h:
                yc = 0.5 * (r["y0"] + r["y1"])
                t.SetStart(V(r["x0"] + h / 2, yc)); t.SetEnd(V(r["x1"] - h / 2, yc)); t.SetWidth(P.FromMM(h))
            else:
                xc = 0.5 * (r["x0"] + r["x1"])
                t.SetStart(V(xc, r["y0"] + w / 2)); t.SetEnd(V(xc, r["y1"] - w / 2)); t.SetWidth(P.FromMM(w))
            t.SetLayer(LAYER[r["layer"]])
            t.SetNet(nets[r["net"]])
            b.Add(t)
            # track ends are round: add a filled rectangle under the ends so the copper
            # outline stays exactly rectangular (identical to the solver geometry)
        z = None
        if min(w, h) >= 0.3:
            # copper pour: solid pad connection, 0.1 mm clearance, no thermal reliefs,
            # min width 0.1 mm, keep islands (all copper is intentional)
            z = P.ZONE(b)
            z.SetLayer(LAYER[r["layer"]])
            z.SetNet(nets[r["net"]])
            ol = z.Outline()
            ol.NewOutline()
            for (x, y) in ((r["x0"], r["y0"]), (r["x1"], r["y0"]), (r["x1"], r["y1"]), (r["x0"], r["y1"])):
                ol.Append(P.FromMM(OX + x), P.FromMM(OY - y))
            z.SetAssignedPriority(len(zones))   # distinct priorities (overlapping same-net pours)
            z.SetLocalClearance(P.FromMM(0.1))
            z.SetPadConnection(P.ZONE_CONNECTION_FULL)
            z.SetMinThickness(P.FromMM(0.1))
            try:
                z.SetIslandRemovalMode(P.ISLAND_REMOVAL_MODE_NEVER)
            except AttributeError:
                pass
            b.Add(z)
            zones.append(z)
    for v in d["vias"]:
        via = P.PCB_VIA(b)
        via.SetPosition(V(v["x"], v["y"]))
        via.SetDrill(P.FromMM(v["drill"]))
        via.SetWidth(P.FromMM(v["dia"]))
        via.SetNet(nets[v["net"]])
        b.Add(via)
    for f in d["fets"]:
        fp = P.FOOTPRINT(b)
        fp.SetReference(f["ref"])
        fp.SetValue(f["value"])
        b.Add(fp)
        for pd in f["pads"]:
            pad = P.PAD(fp)
            try:
                pad.SetShape(P.F_Cu, P.PAD_SHAPE_ROUNDRECT)
            except TypeError:
                pad.SetShape(P.PAD_SHAPE_ROUNDRECT)
            pad.SetAttribute(P.PAD_ATTRIB_SMD)
            pad.SetLayerSet(pad.SMDMask())
            w, h = pd["x1"] - pd["x0"], pd["y1"] - pd["y0"]
            pad.SetSize(P.F_Cu, P.VECTOR2I(P.FromMM(w), P.FromMM(h)))
            pad.SetPosition(V(0.5 * (pd["x0"] + pd["x1"]), 0.5 * (pd["y0"] + pd["y1"])))
            pad.SetNumber(pd["num"])
            pad.SetNet(nets[pd["net"]])
            fp.Add(pad)
    for c in d["caps"]:
        fp = P.FOOTPRINT(b)
        fp.SetReference(c["ref"])
        fp.SetValue(c["value"])
        fp.SetPosition(V(c["x"], c["y"]))   # flip pivot = part centre
        b.Add(fp)
        # 0402: pads 0.5 x 0.5 mm, centres +-0.5 mm along the cap axis (IPC nominal)
        for num, sgn in (("1", -1), ("2", +1)):
            pad = P.PAD(fp)
            pad.SetAttribute(P.PAD_ATTRIB_SMD)
            pad.SetLayerSet(pad.SMDMask())
            pad.SetSize(P.F_Cu, P.VECTOR2I(P.FromMM(0.5), P.FromMM(0.5)))
            if c["rot"] == 0:
                pos = (c["x"] + sgn * 0.5, c["y"])
            else:
                pos = (c["x"], c["y"] + sgn * 0.5)
            # net = net of the copper rectangle (same layer) under the pad centre
            hits = [r["net"] for r in d["rects"] if r["layer"] == c["layer"]
                    and r["x0"] <= pos[0] <= r["x1"] and r["y0"] <= pos[1] <= r["y1"]]
            assert len(set(hits)) == 1, f"cap pad at {pos} lands on {hits}"
            net = hits[0]
            pad.SetPosition(V(*pos))
            pad.SetNumber(num)
            pad.SetNet(nets[net])
            fp.Add(pad)
        if c["layer"] == "B.Cu":
            try:
                fp.Flip(fp.GetPosition(), P.FLIP_DIRECTION_TOP_BOTTOM)
            except (AttributeError, TypeError):
                fp.Flip(fp.GetPosition(), False)
    for fp in b.GetFootprints():          # small, unobtrusive reference text
        for t in (fp.Reference(), fp.Value()):
            t.SetTextSize(P.VECTOR2I(P.FromMM(0.3), P.FromMM(0.3)))
            t.SetTextThickness(P.FromMM(0.05))
            t.SetVisible(False)
    o = d["outline"]
    for (x0, y0, x1, y1) in ((o["x0"], o["y0"], o["x1"], o["y0"]), (o["x1"], o["y0"], o["x1"], o["y1"]),
                             (o["x1"], o["y1"], o["x0"], o["y1"]), (o["x0"], o["y1"], o["x0"], o["y0"])):
        s = P.PCB_SHAPE(b)
        s.SetShape(P.SHAPE_T_SEGMENT)
        s.SetStart(V(x0, y0))
        s.SetEnd(V(x1, y1))
        s.SetLayer(P.Edge_Cuts)
        s.SetWidth(P.FromMM(0.1))
        b.Add(s)
    if zones:
        filler = P.ZONE_FILLER(b)
        filler.Fill(b.Zones())
    P.SaveBoard(dst, b)
    write_stackup(dst, st)
    write_project(dst)


def write_project(board_path):
    """Minimal .kicad_pro with fab-class rules (JLCPCB 4-layer standard capability:
    0.09-0.1 mm min. clearance, 0.3 mm min. drill; check the fab's current page).
    Solder-mask web is set to 0: the EPC2015C land pattern is solder-mask-defined
    at 0.4 mm pitch (datasheet), which a fab handles with a mask-defined opening."""
    import json as _j
    import os
    pro = {
        "board": {"design_settings": {
            "rules": {"min_clearance": 0.1, "min_track_width": 0.1, "min_via_diameter": 0.5,
                      "min_through_hole_diameter": 0.3, "min_copper_edge_clearance": 0.3,
                      "solder_mask_to_copper_clearance": 0.0},
            "rule_severities": {"solder_mask_bridge": "ignore"}}},
        "net_settings": {"classes": [{"name": "Default", "clearance": 0.1, "track_width": 0.15,
                                      "via_diameter": 0.6, "via_drill": 0.3}]},
        "meta": {"filename": os.path.basename(board_path).replace(".kicad_pcb", ".kicad_pro"), "version": 1},
    }
    open(board_path.replace(".kicad_pcb", ".kicad_pro"), "w").write(_j.dumps(pro, indent=2))


def write_stackup(path, st):
    """Insert the physical stackup: L1-(d1)-L2-(core)-L3-(d1)-L4, symmetric."""
    cu, d1, T = st["cu"], st["d1"], st["T"]
    core = T - 2 * (d1 + cu) - cu          # T is L1-to-L4 copper centre distance
    def cl(name, t):
        return f'      (layer "{name}" (type "copper") (thickness {t:.4f}))'
    def dl(i, kind, t):
        return (f'      (layer "dielectric {i}" (type "{kind}") (thickness {t:.4f}) (material "FR4") '
                f'(epsilon_r 4.4) (loss_tangent 0.02))')
    lines = ['    (stackup',
             '      (layer "F.SilkS" (type "Top Silk Screen"))',
             '      (layer "F.Paste" (type "Top Solder Paste"))',
             '      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))',
             cl("F.Cu", cu), dl(1, "prepreg", d1), cl("In1.Cu", cu), dl(2, "core", core),
             cl("In2.Cu", cu), dl(3, "prepreg", d1), cl("B.Cu", cu),
             '      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))',
             '      (layer "B.Paste" (type "Bottom Solder Paste"))',
             '      (layer "B.SilkS" (type "Bottom Silk Screen"))',
             '      (copper_finish "ENIG")',
             '      (dielectric_constraints yes)',
             '    )']
    txt = open(path, encoding="utf-8").read()
    txt = re.sub(r"\(stackup.*?\n\s*\)\n", "", txt, flags=re.S)  # drop a default stackup if present
    txt = txt.replace("(setup", "(setup\n" + "\n".join(lines), 1)
    open(path, "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
