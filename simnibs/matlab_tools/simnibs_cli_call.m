function call = simnibs_cli_call(command)
 % Creates a call to a a SimNBIS CLI command
 call = [SIMNIBSPYTHON ' ' '"' SIMNIBSDIR filesep 'cli' filesep command '.py' '"'];
end