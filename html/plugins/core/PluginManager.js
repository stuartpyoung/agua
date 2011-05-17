
/**
 * OBJECT:  PLUGIN MANAGER, Version 0.1
 * PURPOSE: MANAGE EXTERNAL PLUGINS ON TOP OF dojo FRAMEWORK
 * LICENCE: Copyright (c) 2008 Stuart Young youngstuart@hotmail.com
 *          This code is freely distributable under the terms of an MIT-style license.
 * 
 * METHODS:
 *
 *    PluginManager._installedPlugins(String name)       // INSTALLED PLUGINS
 *    PluginManager.pluginsList()                       // LIST OF PLUGINS (ALT: GET FROM plugins.json FILE)
 *    PluginManager.loadPlugins(String name)            // LOAD PLUGINS FROM LIST
 *    PluginManager.getPluginJson(String name)             // RETURN RETRIEVED JSON FOR PLUGIN
 *    PluginManager.getInstalledPlugins(String name)    // RETURN INSTALLED PLUGINS

 *    PluginManager.getVersion(String name)             // RETURN PLUGIN VERSION
 *    PluginManager.getDescription(String name)         // RETURN PLUGIN DESCRIPTION
 *    PluginManager.getPluginUrl(String name)           // RETURN PLUGIN URL
 *    PluginManager.setInstalled(String name)           // SET PLUGIN INSTALLED Boolean
 *    PluginManager.isInstalled(String name)            // GET PLUGIN INSTALLED Boolean
 *    PluginManager.getInfo(String name)                // GET INFO FOR A PLUGIN
 *
 *        Boolean isInstalled
 *        String  version
 *        String  description
 *        String  pluginPage   URL to download the plugin
 *
 *    DO THESE LATER IF NECESSARY
 *    PluginManager.unInstall(String name)              // REMOVE PLUGIN FROM plugins.json LIST
 *    PluginManager.install(String name)                // ADD PLUGIN TO plugins.json LIST
 *
 * CHANGELOG:
 * Sat 18th October 2008: Version 0.1
 *   load plugins
 *   added license
 * 
 * The comments section is optional, but please leave the copyright
/*--------------------------------------------------------------------*/

dojo.provide("plugins.core.PluginManager");

dojo.require("plugins.core.Plugin");

// OBJECT:  PluginManager
// PURPOSE: LOAD ALL PLUGINS

dojo.declare( "plugins.core.PluginManager", null,
{
	////}

// HASH OF INSTALLED PLUGINS
_installedPlugins : {},

plugins : [],

parentWidget : null,

// LIST OF PLUGINS TO INSTALL
pluginsList :
[
	"plugins.project.Controller"
	,
	"plugins.workflow.Controller"
	,
	"plugins.report.Controller"
	,
	"plugins.view.Controller"
	,
	"plugins.admin.Controller"
],

constructor : function(args) {
// CONSTRUCTOR	

	// SET INPUT PLUGINS LIST IF PROVIDED
	if ( args.pluginsList != null && args.pluginsList )
	{
		this.pluginsList = args.pluginsList;
	}


	// SET PARENT WIDGET IF PROVIDED
	if ( args.parentWidget != null && args.parentWidget )
	{
		this.parentWidget = args.parentWidget;
	}


	// LOAD CORE

	// GET CORE JSON
	var coreJson = this.getPluginJson("plugins.core.Controls");

	// CREATE CORE PLUGIN OBJECT
	var plugin = new plugins.core.Plugin(coreJson);

	// ADD PLUGIN NAME TO INSTALLED PLUGINS HASH
	// ALONG WITH installed STATUS
	this._installedPlugins["plugins.core.Controls"] = plugin;


	//var controls = new plugins.core.Controls( { id: 1, target: 'controls' } );


	// LOAD PLUGINS
	// (ALTERNATELY, GET this.plugins FROM A JSON FILE...)
	this.loadPlugins();
},


loadPlugins : function ()   {
// returns Array of plugin names 


	for ( var i = 0; i < this.pluginsList.length; i++ )
	{
		var pluginName = this.pluginsList[i];

		// SHOW LOADING MESSAGE
//			this.parentWidget.showPreloadMessage("Loading " + pluginName);

		// GET JSON
		var pluginJson = this.getPluginJson(pluginName);

		// CHECK DEPENDENCIES
		if ( this.checkDependencies(pluginJson.dependencies) )
		{
			// INSTALL THE PLUGIN AND AT THE SAME TIME
			// INSTANTIATE A PLUGIN OBJECT FOR IT
			var newPlugin = new plugins.core.Plugin(pluginJson);

			if ( newPlugin != null )
			{
				this.plugins.push(newPlugin);

				// ADD PLUGIN NAME TO INSTALLED PLUGINS HASH
				// ALONG WITH installed STATUS
				this._installedPlugins[pluginJson.name] = newPlugin;

				var moduleName = pluginJson.name.match(/^plugins\.(.+)\./)[1];




				// THIS DOESN'T WORK
				//Agua.controllers[moduleName] = new (eval (pluginJson.name)());

				// INSTANTIATE WIDGET
				Agua.controllers[moduleName] = eval ("new " + pluginJson.name + "()");


				//// DEBUG
				//this.workflowController.createTab();
			}
			else
			{
			}
		}            
	}
},


//getPluginJson : function (pluginName)
//{
//    // GET PATH FOR PLUGINS
//    // E.G.: ../dojo-1.5.0/dojo/
//    var dojoPath = dojo.baseUrl;
//    var jsonPath = dojoPath + "../../";
//
//    // LOAD PLUGIN info.json FILE
//    var infoJson = pluginName.match( /^(.+\.)[^\.]+$/ )[1];
//    infoJson = infoJson.replace( /\./g, "/");
//    infoJson = jsonPath + infoJson + "info.json";    
//    
//    // FETCH workflowTab JSON
//    var json = '';
//    dojo.xhrGet(
//        {
//            url: infoJson,
//            handleAs: "json-comment-optional",
//            sync: true,
//            handle: function(response){
//                json = response;
//            }
//        }
//    );
//
//    return json;  
//},



getPluginJson : function (pluginName) {

	// GET PATH FOR PLUGINS
	// E.G.: ../dojo-1.5.0/dojo/
	var dojoPath = dojo.baseUrl;
	var jsonPath = dojoPath + "../../";

	// LOAD PLUGIN info.json FILE
	//var infoJson = pluginName.match( /^(.+\.)[^\.]+$/ )[1];
	var infoJson = pluginName.match( /^plugins\.(.+\.)[^\.]+$/ )[1];
	infoJson = infoJson.replace( /\./g, "/");
	//var infoJsonUrl = jsonPath + infoJson + "info.json";    
	var infoJsonUrl = dojo.moduleUrl("plugins", infoJson + "info.json");

	// FETCH workflowTab JSON
	var json = '';
	dojo.xhrGet(
		{
			url: infoJsonUrl,
			handleAs: "json-comment-optional",
			sync: true,
			handle: function(response){
				json = response;
			}
		}
	);

	return json;  
},


checkDependencies : function (dependencies) {
	// CHECK DEPENDENCIES ARE ALREADY LOADED AND CORRECT VERSION


	//console.dir(this._installedPlugins);



	if ( ! dependencies )	{	return 1;	}



	for ( var i = 0; i < dependencies.length; i++ )
	{
		var requiredName = dependencies[i].name;
		var requiredVersion = dependencies[i].version;

		// CHECK DEPENDENCY CLASS IS LOADED
		if ( requiredName )
		{
			var dependency = this._installedPlugins[requiredName];

			// CHECK VERSION IS MINIMUM OR GREATER
			if ( dependency.version >= requiredVersion  )    
			{        
				// CHECK THAT THE DEPENDENCY ACTUALLY INSTALLED OKAY
				if ( ! dependency.installed )
				{
					return 0;
				}
				else
				{
				}
			}
			else
			{
				return 0;
			}
		}
		else
		{
			return 0;
		}
	}

	return 1;        
},


// RETURN HASH OF INSTALLED PLUGINS
getInstalledPlugins : function ()   {
	return this._installedPlugins;
}

});

// end of PluginManager

