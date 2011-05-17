dojo.provide( "plugins.core.Plugin");

/**
 * PLUGIN FRAMEWORK, Version 0.1
 * Copyright (c) 2008 Stuart Young youngstuart@hotmail.com
 * This code is freely distributable under the terms of an MIT-style license.
 * 
 *  This code provides the following functions
 *
 *    Plugin.isInstalled(String name)           // RETURN Boolean PLUGIN INSTALL STATUS
 *    Plugin.getVersion(String name)            // RETURN PLUGIN VERSION
 *    Plugin.getDescription(String name)        // RETURN PLUGIN DESCRIPTION
 *    Plugin.getPluginPage(String name)         // RETURN PLUGIN URL
 *    Plugin.getInfo(String name)               // RETURN PLUGIN INFO (NAME, VERSION, DESCRIPTION, IS INSTALLED)
 *
 *        Boolean isInstalled
 *        String  version
 *        String  description
 *        String  pluginPage   URL to download the plugin
 *
 * CHANGELOG:
 * Sat 18th October 2008: Version 0.1
 *   load plugins
 *   load plugins
 *   added license
 * 
 * you may remove the comments section, but please leave the copyright
/*--------------------------------------------------------------------*/



// OBJECT:  Plugin
// PURPOSE: ATTEMPT TO LOAD A PLUGIN USING dojo.require AND STORE
//			WHETHER THE LOAD WAS SUCCESSFUL OR NOT AS installed=BOOLEAN

dojo.declare( "plugins.core.Plugin", null,
{
    name : '',

    version: 0.01,

    description : '',

    url: '',

    dependencies : [],    

    installed : false,


	// CONSTRUCTOR	
	constructor : function(json) {


        if ( ! json.name )
        {
            return;
        }
        this.name = json.name;

        // SET DEPENDENCIES
        this.version = json.version;
        this.description = json.description;
        this.url = json.url;  
        this.dependencies = json.dependencies;        

        var success = this.loadPlugin();
	},

    loadPlugin : function ()
    {
        // LOAD PLUGINS

		if ( this.name != null && this.name && this.name != '' )
		{
			//try {
			//	
			//	var Class = this.name;

				// NB: DON'T USE dojo.require TO AVOID THIS ERROR
				// error loading uri: ./../../../builds/004/dojo/plugins/core/Plugin.js, exception:
				//  TypeError: Cannot call method "split" of undefined


				if ( dojo["require"](this.name) )
				{
					// setInstalled IF LOAD SUCCEEDED
					this.installed = true;
					return 1;
				}
				else
				{
					return 0;
				}
			//}
			//catch (e) {};
		}
    }

	,


    setInstalled : function ()
    {
        if ( this.installed != false && installed != true )
        {
            return 0;
        }

        this.installed = true;
        return 1;
    },

    getInstalled : function ()
    {
        return this.installed;    
    },

    getVersion : function ()
    {
        return this.version;
    },

    getDescription : function ()
    {
        return this.description;
    },

    getPluginUrl : function ()
    {
        return this.pluginUrl;
    },

    getInfo : function ()
    {
        var info = '';
        info += 'Status: ';
        info += this.getInstalled();
        info += '\n';
        info += 'Version: ';
        info += this.version();
        info += '\n';
        info += 'Description: ';
        info += this.description();
        info += '\n';
        info += 'Plugin Url: ';
        info += this.pluginUrl();
        info += '\n';

        return info;
    }

});    
    