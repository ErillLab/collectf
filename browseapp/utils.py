def id2dbxref(id):
    return hex(id)[2:].zfill(5)

def dbxref2id(dbxref):
    return int(dbxref, 16)
