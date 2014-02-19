COLLECTF_TEMPLATE_SBT = '''
Submit-block ::= {
  contact {
    contact {
      name name {
        last "Erill",
        first "Ivan"
      },
      affil std {
        affil "University of Maryland Baltimore County",
        div "Biological Sciences",
        city "Baltimore",
        sub "Maryland",
        country "United States",
        street "1000 Hilltop Circle",
        email "erill@umbc.edu",
        fax "001-410-455-3875",
        phone "001-410-455-2470",
        postal-code "21250"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Erill",
            first "Ivan",
            initials "I.",
            suffix ""
          }
        }
      },
      affil std {
        affil "University of Maryland Baltimore County",
        div "Biological Sciences",
        city "Baltimore",
        sub "Maryland",
        country "United States",
        street "1000 Hilltop Circle",
        postal-code "21250"
      }
    }
  },
  subtype new
}

Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Erill",
              first "Ivan",
              initials "I.",
              suffix ""
            }
          }
        },
        affil std {
          affil "University of Maryland Baltimore County",
          div "Biological Sciences",
          city "Baltimore",
          sub "Maryland",
          country "United States",
          street "1000 Hilltop Circle",
          postal-code "21250"
        }
      },
      title "CollecTF database submission"
    }
  }
}

'''
