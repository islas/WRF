import sane
import itertools

@sane.register
def add_build_for_envs( orch ):
  envs = [ "gnu", "intel-classic", "intel-oneapi", "nvhpc" ]
  sm_opt = [ "ON", "OFF" ]
  dm_opt = [ "ON", "OFF" ]
  build_types = [ "Release", "Debug" ]
  configurations = { "ARW" : [ "EM_REAL", "EM_FIRE", "EM_B_WAVE" ] }
  for core, env, build_type, sm, dm in itertools.product( configurations, envs, build_types, sm_opt, dm_opt ):
    orch.log( f"Creating builds for permutation core: {core} env: {env} build_type: {build_type} sm: {sm} dm: {dm}" ) 
    for case in configurations[core]:
      sm_desc = "_sm" if sm == "ON" else ""
      dm_desc = "_dm" if dm == "ON" else ""
      id = f"{core}_{case}_{env}_{build_type}{sm_desc}{dm_desc}"

      action = sane.Action( f"build_{id}" )
      action.config["command"]   = ".workflow/scripts/buildCMake.sh"
      args = []

      build_dir = "_${{ id }}"
      install_dir = f"install_{id}"

      config_cmd = [ "-p ${{ environment }}", "-d", build_dir, "-i", install_dir ]
      config_cmd.extend( [ "-x -- -DWRF_NESTING=BASIC", f"-DCMAKE_BUILD_TYPE={build_type}" ] )
      config_cmd.extend( [ f"-DWRF_CORE={core}", f"-DWRF_CASE={case}" ] )
      config_cmd.extend( [ f"-DUSE_MPI={dm}", f"-DUSE_OPENMP={sm}" ] )
      build_cmd = [ build_dir, "-j ${{ resources.cpus }}" ]
      clean_cmd = [ "-d", build_dir, "-i", install_dir ]

      args.extend( [ "-c", " ".join( config_cmd ) ] )
      args.extend( [ "-b", " ".join( build_cmd ) ] )
      args.extend( [ "-r", " ".join( clean_cmd ) ] )

      action.config["arguments"] = args
      action.add_resource_requirements( { "cpus" : 8 } )
      action.environment = env
      orch.add_action( action )
