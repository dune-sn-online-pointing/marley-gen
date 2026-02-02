// Example MARLEY job configuration file
// Steven Gardiner <gardiner@fnal.gov>
// Revised 11 May 2020 for version 1.2.0
//
// INTRODUCTION
//
// The MARLEY command-line executable is configured using a JSON-like file
// format. While the file format is quite similar to standard JSON (see
// http://www.json.org/ for a full description), MARLEY job configuration files
// differ from JSON files in the following ways:
//
//   - Single-word keys (no whitespace) may be given without surrounding
//     double quotes
//
//   - C++-style comments // and /* */ are allowed anywhere in the file.
//
//   - A trailing comma is allowed at the end of JSON objects and arrays
//
//  The JSON parser included with MARLEY will print an error message
//  if parsing the job configuration file fails. This message will
//  provide some (hopefully useful) guidance in troubleshooting
//  formatting mistakes. To simplify writing their own configuration files,
//  users are encouraged to copy the examples (particularly
//  "examples/COPY_ME.js") and modify them to suit their needs.
//
//  A file extension of ".js" is recommended for MARLEY configuration files
//  because typical syntax highlighting settings for JavaScript work well with
//  the configuration file format.
//
{ // An opening curly brace begins the configuration file content

    // RANDOM NUMBER SEED (optional)
    //
    // The "seed" key provides a 64-bit unsigned integer (between 0
    // and 2^64 - 1, inclusive) that will be used to seed MARLEY's random
    // number generator.
    //
    // If this key is omitted, MARLEY will use the system time since the
    // Unix epoch as its random number seed.
    seed: 123456,
  
    // INCIDENT NEUTRINO DIRECTION (optional)
    //
    // The "direction" JSON object stores a 3-vector that represents the
    // direction of the incident neutrinos. Note that this vector does not need
    // to be normalized to unity, but at least one element *must* be nonzero.
    //
    // The "x", "y", and "z" keys give the Cartesian components of the direction
    // 3-vector. If the "direction" JSON object is omitted, then x = 0.0, y =
    // 0.0, z = 1.0 will be assumed.
    //
    // An isotropic neutrino direction may be randomly sampled for each event by
    // using the following configuration:
    //
    // direction: "isotropic",
    //
    // In this example, incident neutrinos travel in the +z direction.
    direction: { x: 0.0, y: 0.0, z: 1.0 },
  
    // TARGET SPECIFICATION (optional)
    //
    // The nuclidic composition of the material illuminated by the incident
    // neutrinos in a MARLEY simulation may be specified using the "target"
    // JSON object. This object defines two arrays, which must have equal sizes.
    // The "nuclides" array contains one or more nuclear PDG codes with one entry
    // per distinct nuclide present in the target material. The "atom_fractions"
    // array, as its name suggests, contains the corresponding atom fractions.
    // Elements of the "atom_fractions" array will automatically be normalized
    // to sum to unity if this is not done already by the user. Negative
    // elements in this array will trigger an error message from MARLEY.
    //
    // If the "target" object is omitted from the job configuration file,
    // then every nuclide that appears in the initial state of at least one
    // configured reaction (see the REACTION INPUT FILES(S) section below)
    // will be assumed to be present with equal abundance.
    target: {
      nuclides: [ 1000180400 ],
      atom_fractions: [ 1.0 ],
    },
  
    // REACTION INPUT FILE(S) (required)
    //
    // For simulating neutrino-nucleus reactions, MARLEY relies on tables of
    // precomputed nuclear matrix elements to compute neutrino reaction cross
    // sections. The nuclear matrix elements to use in any particular simulation
    // job are specified using a JSON array stored under the "reactions" key.
    // Each entry in the array is the name of a data file that contains a set of
    // nuclear matrix elements to use when calculating cross sections for a
    // particular scattering mode on a particular target nuclide. Note that the
    // file names must be simple strings; MARLEY configuration files currently do
    // not support the use of environment variables, bash globs, etc.
  
    // MARLEY also uses an input file to configure neutrino-electron
    // elastic scattering (ES) reactions. In this case, the file provides a list
    // of the nuclides for which the ES process should be simulated.
  
    // Although multiple reaction input files may be available for a particular
    // process, only a single configuration is allowed for any particular
    // combination of target nucleus and scattering mode. Conflicting
    // configurations will be ignored in favor of the one that appears in the
    // first relevant reaction input file listed in the "reactions" array.
    //
    // In MARLEY v1.2.0, the main process available to be simulated is
    // charged-current scattering of electron neutrinos on 40Ar.
    // There are three available reaction input files (stored in the
    // folder data/react/) for this process
    //
    //   - ve40ArCC_Bhattacharya2009.react
    //   - ve40ArCC_Bhattacharya1998.react
    //   - ve40ArCC_Liu1998.react
    //
    // See each of these files for details about their respective nuclear matrix
    // element evaluations.
    //
    // Two other reaction input files are currently included in the official
    // MARLEY source code distribution:
    //
    //   - ES.react: Enables simulation of neutrino-electron elastic scattering
    //               on a 40Ar atomic target
    //
    //   - CEvNS40Ar.react: Enables simulation of coherent elastic
    //                      neutrino-nucleus scattering on a 40Ar target. The
    //                      q^2 dependence of the nuclear form factor is
    //                      neglected.
    //
    // The full path to the file does not need to be given in each element of the
    // "reactions" JSON array. After searching in the working directory, MARLEY's
    // default behavior is to search for reaction input files in the directories
    // ${MARLEY}/data, ${MARLEY}/data/react/, and ${MARLEY}/data/structure/,
    // where ${MARLEY} is the value of the MARLEY environment variable (typically
    // set by sourcing the setup_marley.sh bash script). The list of search
    // directories beyond the working directory can be changed from the default
    // by setting the MARLEY_SEARCH_PATH environment variable to a ':'-separated
    // list of directories.
    //
    // Unless a particular reaction channel is represented by a data file given
    // in the "reactions" JSON array, it will not be included in the MARLEY
    // simulation.
    //
    reactions: [ "ES.react" ],
  
    // NEUTRINO SOURCE SPECIFICATION (required)
    //
    // The "source" JSON object describes the incident neutrino spectrum. The
    // input spectrum should not be weighted by a reaction cross section because
    // MARLEY will weight the incident spectrum using its own cross section
    // model during the simulation. The flux units used to specify user-defined
    // spectra will be ignored by MARLEY (the source spectrum will be
    // renormalized to unity internally to produce a probability density
    // function), but all neutrino energies should be given in MeV.
    //
    // The "type" key owned by the source specification describes the
    // incident spectrum and determines the other parameters that must be
    // specified. The currently allowed values for the "type" key are
    // given in the following table (with synonyms separated by commas)
    //
    //   Spectrum description                   Allowed "type" values
    //   --------------------                   ---------------------
    //
    //   Fermi-Dirac                            "fd", "fermi-dirac",
    //                                          "fermi_dirac"
    //
    //
    //   "Beta-fit"                             "bf", "beta", "beta-fit"
    //   (see, e.g.,
    //   http://arxiv.org/abs/1511.00806)
    //
    //   Monoenergetic                          "mono", "monoenergetic"
    //
    //   Muon decay-at-rest                     "dar", "decay-at-rest"
    //   (ve and vu only)
    //
    //   User-defined histogram                 "hist", "histogram"
    //
    //   User-defined probability               "grid"
    //   density function evaluated via
    //   interpolation on a set of grid
    //   points
    //
    //   ROOT TH1                               "th1"
    //
    //   ROOT TGraph                            "tgraph"
    //
    // All of the source types also require the use of the "neutrino" key,
    // which specifies the species of neutrino emitted by the source. Valid
    // values for the "neutrino" key are neutrino PDG codes (±12, ±14, ±16)
    // and the strings "ve", "vebar", "vu", "vubar", "vt", and "vtbar".
    //
    // Although MARLEY's algorithm for sampling from arbitrary spectra is
    // reasonably robust, user-defined spectra with widely-separated tight
    // peaks or other unusual shapes may experience problems. Users should
    // perform some simple verification simulations when using a "hist",
    // "grid", "th1", or "tgraph" neutrino source. For assistance with
    // this testing or to report bugs, please contact the MARLEY developers
    // (support@marleygen.org).
    //
    // Examples of each of the allowed source specifications are shown below
    //
    // FERMI-DIRAC
    //
    //  source: {
    //    type: "fermi-dirac",
    //    neutrino: "ve",
    //    Emin: 0,           // Minimum neutrino energy (MeV)
    //    Emax: 60,          // Maximum neutrino energy (MeV)
    //    temperature: 3.5,  // Temperature (MeV)
    //    eta: 4             // Pinching parameter (dimensionless, default 0)
    //  },
    //
    // "BETA FIT"
    //
    //  source: {
    //    type: "beta-fit",
    //    neutrino: "ve",
    //    Emin: 0,           // Minimum neutrino energy (MeV)
    //    Emax: 60,          // Maximum neutrino energy (MeV)
    //    Emean: 15,         // Mean neutrino energy (MeV)
    //    beta: 3.0,         // Pinching parameter (dimensionless, default 4.5)
    //  },
    //
    //  MONOENERGETIC
    //
    //  source: {
    //    type: "monoenergetic",
    //    neutrino: "ve",
    //    energy: 10,        // Neutrino energy (MeV)
    //  },
    //
    //  MUON DECAY-AT-REST
    //
    //  source: {
    //    type: "decay-at-rest",
    //    neutrino: "ve",
    //  },
    //
    //  HISTOGRAM
    //
    //  source: {
    //    type: "histogram",
    //    neutrino: "ve",
    //    E_bin_lefts: [ 10., 20., 30. ],   // Low edges of energy bins (MeV)
    //    weights: [ 0.2, 0.5, 0.3 ],       // Bin weights (dimensionless)
    //    Emax: 40.,                        // Upper edge of the final bin (MeV)
    //  },
    //
    //  Within a histogram bin, energies are sampled uniformly on the
    //  half-open interval [ Ebin_left, Ebin_right ).
    //
    //  GRID
    //
    //  source: {
    //    type: "grid",
    //    neutrino: "ve",
    //    energies: [ 10., 15., 20. ],   // Energy grid points (MeV)
    //
    //    prob_densities: [ 0., 1., 0. ],  // Probability densities
    //                                     // (dimensionless, do not need to be
    //                                     //  normalized to unity by the user)
    //
    //    rule: "linlin",                // Interpolation rule ("linlin" default)
    //  },
    //
    //  The allowed values of the "rule" key are "linlin" (linear-linear
    //  interpolation), "loglog" (log-log interpolation), "linlog" (linear
    //  in energy, logarithmic in probability density), and "loglin"
    //  (logarithmic in energy, linear in probability density).
    //
    //
    //  TH1
    //
    //  source: {
    //    type: "th1",
    //    neutrino: "ve",
    //    tfile: "my_root_file.root",  // Name of the ROOT file containing
    //                                 // the TH1 object
    //
    //    namecycle: "MyHist",         // Name under which the TH1 object
    //                                 // appears in the file (used to
    //                                 // retrieve it)
    //  },
    //
    //  TGRAPH
    //
    //  source: {
    //    type: "tgraph",
    //    neutrino: "ve",
    //    tfile: "my_root_file.root",  // Name of the ROOT file containing
    //                                 // the TGraph object
    //
    //    namecycle: "MyGraph",        // Name of the TGraph object (used to
    //                                 // retrieve it from the ROOT file)
    //  },
    //
    //
    //
  
      //  HISTOGRAM
  //
    source: {
       type: "histogram",
       neutrino: "ve",
       E_bin_lefts: [ 2.0 ],   // Low edges of energy bins (MeV)
       weights: [ 1.0 ],       // Bin weights (dimensionless)
       Emax: 70.0,                        // Upper edge of the final bin (MeV)
       weight_flux: false
     },
  
    // EXECUTABLE SETTINGS (optional)
    //
    // The entries within the executable_settings JSON object are used to
    // control the marley command-line executable. They are ignored if
    // the job configuration file is used to initialize MARLEY outside of that
    // context (e.g., within a Geant4 application that links to the MARLEY
    // shared libraries).
    executable_settings: {
  
      // EVENT COUNT (optional)
      // If this key is omitted, a value of 1000 will be assumed.
      events: 1000000,
  
      output: [ { file: "flat_ES.ascii", format: "ascii", mode: "overwrite" } ],
    },
  
  } // A closing curly brace should appear at the end of the file