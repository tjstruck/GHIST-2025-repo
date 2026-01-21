profile_to_username = { 3446506: "@srong", 3473724: "@tjstruck", 3474352: "@RyanGutenkunst", 3481948: "@jopie", 3493671: "@dschride", 3508154: "@austin.t.daigle", 3508270: "@RGutenkunst", 3508906: "@avaughn", 3510655: "@milesdavidroberts", 3511678: "@schien", 3511690: "@resplin5072", 3511775: "@zhaobohu2", 3512651: "@atlashrike", 3514806: "@tsbertino", 3515438: "@noscode", 3515829: "@wenjie.zhu", 3515974: "@jterhorst", 3538708: "@nzuppas", 3545374: "@zoehert", 3547404: "@amhughes", 3547508: "@akwakye", 3547576: "@frayerme13", 3547579: "@ProManna", 3547664: "@as5273", 3547675: "@dragonHunter", 3547678: "@Kadenwinspear", 3547682: "@yamaly", 3547683: "@kiranbajaj", 3547684: "@ekhowell", 3547686: "@celinetree", 3547688: "@ajayior", 3547690: "@paruljohri", 3547691: "@solomonsloat", 3547695: "@KemoTherapy", 3547696: "@camilachavez", 3547697: "@mcatto3", 3547698: "@jacobimarsh", 3547983: "@Alan-Izarraras", 3549066: "@chimarokeonyeaghala", 3549261: "@leann.g.ward", 3550249: "@jannahsaid", 3550636: "@Igelkott", 3550667: "@Ronkyy", 3551012: "@kii925", 3553477: "@Dale_Decena", 3553620: "@alouette", 3553843: "@lzong", 3553877: "@eppleym", 3554268: "@rgollnisch", 3554444: "@wang0207", 3555290: "@currocam", 3556682: "@zmeziere", 3556685: "@Showitt", 3556793: "@Froggolo", 3557522: "@JiatongLiang", 3557947: "@huangdaxian", 3558293: "@nomis-c", 3558940: "@jmurgamoreno", 3559048: "@Mimi33x33x", 3559911: "@tstentella", 3560063: "@arndt", 3561154: "@peterlaurin", 3561587: "@Brendan_A" }

from pathlib import Path
import shutil, glob

for fname in glob.glob("GHIST_submissions_and_truth/final_submissions/*/*"):
    path = Path(fname)
    profile_id = int(path.stem.split("_")[-1])
    if profile_id in profile_to_username:
        username = profile_to_username[profile_id]
        new_fname = fname.replace(f"{profile_id}", f"{username[1:]}")
        print(fname, "->", new_fname)
        shutil.copy2(fname, new_fname)