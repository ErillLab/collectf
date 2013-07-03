def id2dbxref(id):
    return '00' + hex(id)[2:].zfill(5) + '0'

def dbxref2id(dbxref):
    return int(dbxref[2:-1], 16)
