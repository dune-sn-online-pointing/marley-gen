// Use this example job configuration file as a starting point for your own
// files.
{
  seed: 123456, // Random number seed (omit to use time since Unix epoch)

  // Pure 40Ar target
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  // Simulate CC ve scattering on 40Ar
  reactions: [ "ve40ArCC_Bhattacharya2009.react" ],

  // Neutrino source specification
  source: {
    neutrino: "ve",        // The source produces electron neutrinos
    type: "monoenergetic",
    energy: 30.0,          // MeV
  },


  // Incident neutrino direction 3-vector
  direction: { x: 0.0, y: 0.0, z: 1.0 },

  // Settings for marley command-line executable
  executable_settings: {

    // The number of events to generate
    events: 10000,

    // Event output configuration
    output: [ { file: "events30.ascii", format: "ascii", mode: "overwrite" } ],

  },
}
