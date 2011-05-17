dojo.provide("plugins.view.Chromosomes");


// DISPLAY DIFFERENT PAGES TO ALLOW THE view AND ORDINARY
// USERS TO MODIFY THEIR SETTINGS

//// DnD
//dojo.require("dojo.dnd.Source"); // Source & Target
//dojo.require("dojo.dnd.Moveable");
//dojo.require("dojo.dnd.Mover");
//dojo.require("dojo.dnd.move");
//
//// fisheye lite menu animation
//dojo.require("dojox.widget.FisheyeLite");
//
//// comboBox data store
//dojo.require("dojo.data.ItemFileReadStore");
//dojo.require("dijit.form.ComboBox");
//
//// rightPane buttons
//dojo.require("dijit.form.Button");



// STORE FOR PROJECT AND WORKFLOW COMBOS
dojo.require("dojo.data.ItemFileReadStore");

// CUSTOM EDITOR
dojo.require("plugins.report.Editor");


// FILE UPLOAD
dojo.require("plugins.upload.FileUpload");

// RESIZE
dojo.require("dojox.layout.ResizeHandle");

// NOTES EDITOR
dojo.require("dijit.Editor");
dojo.require("dijit.form.DateTextBox");
dojo.require("dijit.form.Textarea");

dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.ValidationTextBox");
dojo.require("dijit.form.NumberTextBox");
dojo.require("dijit.form.CurrencyTextBox");
dojo.require("dojo.currency");
dojo.require("dijit.Dialog");


// FILE MANAGER HAS FILE SELECTORS
dojo.require("plugins.workflow.FileManager");
dojo.require("plugins.workflow.FileSelector");


dojo.require("dojo.data.ItemFileReadStore");
dojo.require("dijit.form.ComboBox");
dojo.require("dijit.Tree");
dojo.require("dijit.layout.AccordionContainer");
dojo.require("dijit.layout.TabContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.FloatingPane");
dojo.require("dojo.fx.easing");
dojo.require("dojox.rpc.Service");
dojo.require("dojo.io.script");


// INHERITED
dojo.require("plugins.core.Common");

// HAS A
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.ExpandoPane");



// WIDGETS IN TEMPLATE
dojo.require("dijit.layout.SplitContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("dojo.parser");


// JBROWSE
//dojo.require("plugins.view.jbrowse.jslib.dojo.jbrowse_dojo");
dojo.require("plugins.view.jbrowse.js.Browser");
dojo.require("plugins.view.jbrowse.js.Util");
dojo.require("plugins.view.jbrowse.js.NCList");
dojo.require("plugins.view.jbrowse.js.Layout");
dojo.require("plugins.view.jbrowse.js.LazyArray");
dojo.require("plugins.view.jbrowse.js.LazyPatricia");
dojo.require("plugins.view.jbrowse.js.Track");
dojo.require("plugins.view.jbrowse.js.SequenceTrack");
dojo.require("plugins.view.jbrowse.js.FeatureTrack");
dojo.require("plugins.view.jbrowse.js.UITracks");
dojo.require("plugins.view.jbrowse.js.ImageTrack");
dojo.require("plugins.view.jbrowse.js.GenomeView");
dojo.require("plugins.view.jbrowse.js.touchJBrowse");

// SPECIES- AND BUILD-SPECIFIC FILES
dojo.require("plugins.view.jbrowse.species.human.hg19.data.refSeqs");
dojo.require("plugins.view.jbrowse.species.human.hg19.data.trackInfo");
//dojo.require("plugins.view.jbrowse.human.data.refSeqs");
//dojo.require("plugins.view.jbrowse.human.data.trackInfo");


dojo.require("plugins.view.jbrowse.js.Browser");

dojo.require("dijit._base.place");

dojo.declare( "plugins.view.Chromosomes",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	////}

// PATH TO WIDGET TEMPLATE
templatePath: dojo.moduleUrl("", "../plugins/view/templates/view.html"),

// PARENT NODE, I.E., TABS NODE
parentWidget : null,

// PROJECT NAME AND WORKFLOW NAME IF AVAILABLE
project : null,
workflow : null,

// onChangeListeners : Array. LIST OF COMBOBOX ONCHANGE LISTENERS
onChangeListeners : [],

// setListeners : Boolean. SET LISTENERS FLAG 
setListeners : false,

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// CSS FILES
cssFiles : ["plugins/view/css/view.css", "plugins/view/css/genome.css"],

// browsers Object. HASH ARRAY OF OPENED BROWSERS
browsers : [],

constructor : function(args) {	

	// NB: TAB HIERARCHY IS AS FOLLOWS:
	//
	//	tabs	
	//		mainTab
	//				leftPane
	//						comboBoxes
	//				rightPane
	//						Browser
	//								Features
	//								GenomeView

	// SET ARGS
	this.parentWidget = args.parentWidget;
	this.project = args.project;
	this.workflow = args.workflow;

	// LOAD CSS
	this.loadCSS();		
},

postMixInProperties: function() {
},

postCreate: function() {

	this.startup();
},


// STARTUP
startup : function () {

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

	// ADD THE PANE TO THE TAB CONTAINER
	this.attachWidget.addChild(this.mainTab);
	this.attachWidget.selectChild(this.mainTab);

	// EXPAND LEFT PANE (WAS CLOSED SO THAT RIGHT PANE WOULD RENDER)
	this.leftPane.toggle();

	// LOAD COMBOS IN SUCCESSION
	this.setProjectCombo();

	// SET COMBO LISTENERS
	setTimeout(
		function (thisObj) { thisObj.setComboListeners(); },
		1000,
		this
	);

	//// LOAD SPARKLINES ON CHROMOSOMES IN MIDDLE PANE
	//this.sparklines();
},

// SHOW SPARK LINES OF EXPRESSION, SNPS, ALIGNMENT COVERAGE, ETC.
// http://localhost/aqwa/0.4/dojo.1.2.2/dojox/charting/tests/test_sparklines.html	
sparklines : function () {

},


// GET FUNCTIONS	
getProject : function () {
	return this.projectCombo.get('value');
},
getWorkflow : function () {
	return this.workflowCombo.get('value');
},
getView : function () {

	return this.viewCombo.get('value');
},
getBuild : function () {

	var buildName;
	var speciesBuild = this.speciesCombo.get('value');
	if ( speciesBuild.match(/^(\S+)\s+\(([^\)]+)\)$/) )
	{
		buildName = speciesBuild.match(/^(\S+)\s+\(([^\)]+)\)$/)[2];
	}

	return buildName;
},
getSpecies : function () {

	var speciesName;
	var speciesBuild = this.speciesCombo.get('value');
	if ( speciesBuild.match(/^(\S+)\s+\(([^\)]+)\)$/) )
	{
		speciesName = speciesBuild.match(/^(\S+)\s+\(([^\)]+)\)$/)[1];
	}

	return speciesName;
},

// SET FUNCTIONS
setProjectCombo : function (projectName, workflowName) {

	// DO COMBO WIDGET SETUP	
	this.inherited(arguments);


	// SET CSS
	this.projectCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.projectCombo.wrapperClass = "view dijitPopup";
	this.projectCombo.itemHeight = 30;		

	if ( projectName == null )	projectName = this.projectCombo.getValue();

	// RESET THE WORKFLOW COMBO
	this.setWorkflowCombo(projectName, workflowName);
},

setWorkflowCombo : function (projectName, workflowName) {

	// DO COMBO WIDGET SETUP	
	this.inherited(arguments);


	// SET CSS
	this.workflowCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.workflowCombo.wrapperClass = "view dijitPopup";
	this.workflowCombo.itemHeight = 30;		

	if ( projectName == null ) projectName = this.projectCombo.getValue();
	if ( workflowName == null ) workflowName = this.workflowCombo.getValue();

	// RESET THE VIEW COMBO
	this.setViewCombo(projectName, workflowName);
},


setViewCombo : function (projectName, workflowName, viewName) {

	// DO COMBO WIDGET SETUP	
	this.inherited(arguments);

	// SET CSS
	this.viewCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.viewCombo.wrapperClass = "view dijitPopup";
	this.viewCombo.itemHeight = 30;		

	// SET DROP TARGET (LOAD MIDDLE PANE, BOTTOM)
	if ( viewName == null ) viewName = this.getView();


	this.setSpeciesCombo(projectName, workflowName, viewName);
},

setSpeciesCombo : function (projectName, workflowName, viewName, speciesName, buildName) {

	// SET CSS
	this.viewCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.viewCombo.wrapperClass = "view dijitPopup";
	this.viewCombo.itemHeight = 30;		

	// SET DROP TARGET (LOAD MIDDLE PANE, BOTTOM)
	if ( viewName == null ) viewName = this.viewCombo.get('value');

	var views = Agua.getViewsByProjectWorkflow(projectName, workflowName);   

	// TAKE THE FIRST VIEW FOR THIS WORKFLOW IF speciesName NOT DEFINED
	if ( views.length != null )
	{
		if ( speciesName == null || ! speciesName )  
		{

			speciesName = views[0].species;
			buildName = views[0].build;

			if ( speciesName == "" ) speciesName = "human";
			if ( buildName == "" ) buildName = "hg19";
		}
	}

	// IF NO VIEWS FOR THIS WORKFLOW,
	// DO NOTHING
	//////////// SELECT FIRST ENTRY IN SPECIES COMBO
	else
	{
		if ( speciesName == null || ! speciesName )  
		{

			//this.speciesCombo.selectedIndex = 2;
		}
	}

	if ( speciesName == null ) speciesName = this.getSpecies();
	if ( buildName == null ) buildName = this.getBuild();


	if ( this.setListeners == false )
	{

		// SET LISTENERS
		setTimeout(function(thisObj) { thisObj.setComboListeners(); }, 2000, this);

		// SET setListeners FLAG
		this.setListeners = true;
	}


	this.loadBrowser(projectName, workflowName, viewName, speciesName, buildName);
},



// LOAD JBROWSE
loadBrowser : function (project, workflow, view, species, build) {
	if ( species == null )	species = this.getSpecies();
	if ( build == null )	build = this.getBuild();


	// SET DEFAULT SPECIES TO HUMAN
	if ( species == null || ! species ) species = "human";

	if ( project == null
		|| workflow == null
		|| view == null
		|| species == null
		|| build == null )
	{
		return;
	}

	var viewObject = Agua.getView(project, workflow, view, species, build);
	var location	=	viewObject.location;
	var trackList	=	viewObject.tracklist;

	// GET UNIQUE ID FOR THIS MENU TO BE USED IN DND SOURCE LATER
	var objectName = "plugins.view.View.jbrowse.Browser";
	var browserId = dijit.getUniqueId(objectName.replace(/\./g,"_"));

	var dataRoot = "plugins/view/jbrowse/" + species + "/";


	//var b = new Browser({
	var b = new plugins.view.jbrowse.js.Browser({
		projectName: project,
		workflowName: workflow,
		viewName: view,
		speciesName: species,
		buildName: build,

		refSeqs: refSeqs,
		trackData: trackInfo,
		baseUrl : "plugins/view/jbrowse/",
		dataRoot : "plugins/view/jbrowse/species/" + species + "/" + build + "/",
		browserRoot : "plugins/view/jbrowse/",
		defaultLocation : "chr1:1..10000000",

		trackList : "Affy_Exon_Probes,Microsatellite,Encode_Common_CNV",

		attachWidget : this.rightPane
		//dataRoot : "plugins/view/jbrowse/users/syoung/Project1/Workflow1/"
	});


	// ADD TO this.browsers ARRAY		
	this.addBrowser(b, project, workflow, view, build, species);		

	//b.showTracks("DNA,gene,mRNA,noncodingRNA");
	//b.navigateTo("chr1:0-20000000");

	if ( location != null && location != '' )
	{
		b.navigateTo(location);
	}

}, // 	loadBrowser 


reloadBrowser : function (project, workflow, view, species, build) {
	if ( species == null )	species = this.getSpecies();
	if ( build == null )	build = this.getBuild();


	if ( project == null
		|| workflow == null
		|| view == null
		|| species == null
		|| build == null )
	{
		return;
	}

	if ( this.getBrowser(project, workflow, view ) )
	{
		// REMOVE EXISTING BROWSER FOR THIS VIEW
		this.removeBrowser(project, workflow, view);
	}

	this.loadBrowser(project, workflow, view, species, build);


}, // 	reloadBrowser 


addBrowser : function (browser, project, workflow, view) {

	var key = project + "*" + workflow + "*" + view;
	this.browsers[key] = browser;
},


getBrowser : function (project, workflow, view) {

	var key = project + "*" + workflow + "*" + view;
	var browser = this.browsers[key];

	return browser;
},


removeBrowser : function (project, workflow, view) {


	var browser = this.getBrowser(project, workflow, view);

	// DESTROY WIDGET
	//browser.destroyRecursive();
	this.rightPane.removeChild(browser.mainTab);

	var key = project + "*" + workflow + "*" + view;
	delete this.browsers[key];
},


deleteView : function () {

	var project = this.getProject();
	var workflow = this.getWorkflow();
	var view = this.getView();



},

setComboListeners : function () {

	// PROJECT COMBO
	var thisObj = this;
	dojo.connect(this.projectCombo, "onchange", dojo.hitch(this, function(projectName) {
			thisObj.setComboListeners(projectName);
	}));

	// WORKFLOW COMBO
	var thisObj = this;
	dojo.connect(this.workflowCombo, "onchange", dojo.hitch(this, function(workflowName) {
			var projectName = thisObj.projectCombo.getValue();
			thisObj.setComboListeners(projectName, workflowName);
	}));

	// VIEW COMBO
	dojo.connect(this.viewCombo, "onchange", dojo.hitch(function(viewName) {

		var projectName = thisObj.projectCombo.getValue();
		var workflowName = thisObj.workflowCombo.getValue();
		thisObj.loadBrowser(projectName, workflowName, viewName);
	}));


	// SPECIES COMBO
	dojo.connect(this.speciesCombo, "onchange", dojo.hitch(function(speciesBuild) {

		var projectName = thisObj.projectCombo.getValue();
		var workflowName = thisObj.workflowCombo.getValue();
		var viewName = thisObj.viewCombo.getValue();

		var speciesName = thisObj.getSpecies();
		var buildName = thisObj.getBuild();

		thisObj.reloadBrowser(projectName, workflowName, viewName, speciesName, buildName);
	}));


	// SET NEW PROJECT LISTENER
	var thisObject = this;
	this.viewCombo._onKeyPress = function(event){

		// summary: handles keyboard events
		var key = event.charOrCode;			
		if ( key == 13 )
		{
			thisObject.workflowCombo._hideResultList();

			var projectName = thisObject.projectCombo.get('value');
			var workflowName = thisObject.workflowCombo.get('value')
			var viewName = thisObject.viewCombo.get('value')

			// STOP PROPAGATION
			event.stopPropagation();

			var isView = Agua.isView(projectName, workflowName, viewName);
			if ( isView == false )
			{
				var viewObject = new Object;
				viewObject.project = projectName;
				viewObject.workflow = workflowName;
				viewObject.name = viewName;

				Agua.addView(viewObject);
				thisObject.setSpeciesCombo(projectName, workflowName, viewName);
			}

			if ( thisObject.viewCombo._popupWidget != null )
			{
				thisObject.viewCombo._showResultList();
			}
		}

		// STOP PROPAGATION
		event.stopPropagation();
	};

}



}); // end of plugins.view.View

