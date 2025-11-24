{
    seed: 123456,
  
    // Neutrino direction: x=1, y=1, z=1 (MARLEY will normalize this)
    direction: { x: 1.0, y: 1.0, z: 1.0 },
  
    target: {
      nuclides: [ 1000180400 ],  // Argon-40
      atom_fractions: [ 1.0 ],
    },
  
    // Elastic scattering only
    reactions: [ "ES.react" ],
  
    // Flat energy spectrum from 2 to 70 MeV
    source: {
       type: "histogram",
       neutrino: "ve",
       E_bin_lefts: [ 2.0 ],
       weights: [ 1.0 ],
       Emax: 70.0,
       weight_flux: false
     },
  
    executable_settings: {
      events: 400,
      output: [ { file: "data/directional_ES.root", format: "root", mode: "overwrite" } ],
    },
}
