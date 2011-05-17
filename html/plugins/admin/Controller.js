dojo.provide("plugins.admin.Controller");

// OBJECT:  plugins.admin.Controller
// PURPOSE: GENERATE AND MANAGE Admin PANES

// INHERITS
//dojo.require("plugins.core.WidgetFramework");

// GLOBAL ADMIN CONTROLLER VARIABLE
var adminController;

// HAS
dojo.require("plugins.admin.Admin");

dojo.declare( "plugins.admin.Controller",
	[ dijit._Widget, dijit._Templated ],
{
	// PANE ID 
	paneId : null,

	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/controller.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// CONSTRUCTOR	
	constructor : function(args) {

		// LOAD CSS FOR BUTTON
		this.loadCSS();

	},


	postCreate : function()
	{

		this.startup();
	},


	startup : function ()
	{

		this.inherited(arguments);

		// ADD ADMIN TAB TO TAB CONTAINER		

		this.adminButton = Agua.addToolbarButton("Admin");


		Agua.toolbar.addChild(this.adminButton);

		// SET ADMIN BUTTON LISTENER
		var listener = dojo.connect(this.adminButton, "onClick", this, "createTab");

		// CREATE TAB
		this.createTab();		

	},


	createTab : function ()
	{

		var admin = new plugins.admin.Admin( { paneId: "adminPane" } );
		this.paneId++;

		// GENERATE DYNAMIC DHTML MENU
		var dynamicMenu = this.createMenu();

		////// BIND THE MENU TO THE DND NODES SO WE OPEN THE
		////// MENU WHEN WE CLICK OVER THE NODES
		////dynamicMenu.bindDomNode( dojo.byId(admin.leftPaneId) );
		////dynamicMenu.bindDomNode( dojo.byId(admin.rightPaneId) );
	},


	// ADD PROGRAMMATIC CONTEXT MENU
	createMenu : function ()
	{
		var dynamicMenu = new dijit.Menu( { id: "admin" + this.paneId + 'dynamicMenuPopup'} );

		// ADD MENU TITLE
		dynamicMenu.addChild(new dijit.MenuItem( { label:"Application Menu", disabled:false} ));
		dynamicMenu.addChild(new dijit.MenuSeparator());

		//// ONE OF FOUR WAYS TO DO MENU CALLBACK WITH ACCESS TO THE MENU ITEM AND THE CURRENT TARGET 	
		// 4. dojo.connect CALL
		//	REQUIRES:
		//		ADDED menu.currentTarget SLOT TO dijit.menu
		var mItem1 = new dijit.MenuItem(
			{
				id: "admin" + this.paneId + "remove",
				label: "Remove",
				disabled: false
			}
		);
		dynamicMenu.addChild(mItem1);
		dojo.connect(mItem1, "onClick", function()
			{
				var parentNode = dynamicMenu.currentTarget.parentNode;
				parentNode.removeChild(dynamicMenu.currentTarget);	
			}
		);

		// SEPARATOR
		dynamicMenu.addChild(new dijit.MenuSeparator());

		//	ADD run MENU ITEM
		var mItem2 = new dijit.MenuItem(
			{
				id: "admin" + this.paneId + "run",
				label: "Run",
				disabled: false
			}
		);
		dynamicMenu.addChild(mItem2);	

		dojo.connect(mItem2, "onClick", function()
			{
				var currentTarget = dynamicMenu.currentTarget; 
				var adminList = currentTarget.parentNode;
			}
		);

		return dynamicMenu;
	},

	loadCSS : function()
	{
		// LOAD CSS
		var cssFiles = [ "plugins/admin/css/controller.css" ];
		for ( var i in cssFiles )
		{
			var cssFile = cssFiles[i];
			var cssNode = document.createElement('link');
			cssNode.type = 'text/css';
			cssNode.rel = 'stylesheet';
			cssNode.href = cssFile;
			cssNode.media = 'screen';
			document.getElementsByTagName("head")[0].appendChild(cssNode);
		}
	}

}); // end of Controller

dojo.addOnLoad(
	function()
	{

		//if ( Agua.loginController == null )
		//{
			//Agua.adminController = new plugins.admin.Controller( { id: 'adminController' } );

			// DEBUG
			//Agua.adminController.createTab();
		//}		
	}
);

