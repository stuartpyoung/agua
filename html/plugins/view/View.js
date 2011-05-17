dojo.provide("plugins.view.View");

/* CLASS SUMMARY: CREATE AND MODIFY AND VIEWS


	TAB HIERARCHY IS AS FOLLOWS:

		tabs	
			mainTab
					leftPane
							comboBoxes
					rightPane
							Browser
									Features
									GenomeView
*/	

// STORE FOR PROJECT AND WORKFLOW COMBOS
dojo.require("dojo.data.ItemFileReadStore");

// CUSTOM EDITOR
//dojo.require("plugins.report.Editor");

// FILE UPLOAD
dojo.require("plugins.upload.FileUpload");

// RESIZE
dojo.require("dojox.layout.ResizeHandle");

// NOTES EDITOR
//dojo.require("dijit.Editor");
//dojo.require("dijit.form.DateTextBox");
//dojo.require("dijit.form.Textarea");

//dojo.require("dijit.form.TextBox");
//dojo.require("dijit.form.ValidationTextBox");
//dojo.require("dijit.form.NumberTextBox");
//dojo.require("dijit.form.CurrencyTextBox");
//dojo.require("dojo.currency");
dojo.require("dijit.Dialog");

// FILE MANAGER HAS FILE SELECTORS
//dojo.require("plugins.workflow.FileManager");
//dojo.require("plugins.workflow.FileSelector");

// INHERITED
dojo.require("plugins.core.Common");

// HAS A
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.ExpandoPane");

// WIDGETS IN TEMPLATE
dojo.require("dijit.layout.SplitContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("dojo.data.ItemFileReadStore");
dojo.require("dijit.form.ComboBox");
dojo.require("dijit.layout.TabContainer");
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.FloatingPane");
dojo.require("dojo.fx.easing");
dojo.require("dojo.parser");
//dojo.require("dijit.layout.ContentPane");
//dojo.require("dijit.Tree");
//dojo.require("dijit.layout.AccordionContainer");
//dojo.require("dojox.rpc.Service");
//dojo.require("dojo.io.script");

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
//dojo.require("plugins.view.jbrowse.species.human.hg19.data.refSeqs");
//dojo.require("plugins.view.jbrowse.species.human.hg19.data.trackInfo");
var refSeqs;
var trackInfo;

dojo.require("dijit._base.place");

dojo.declare( "plugins.view.View",
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
onChangeListeners : new Object,

// setListeners : Boolean. SET LISTENERS FLAG 
setListeners : false,

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// CSS FILES
cssFiles : ["plugins/view/css/view.css", "plugins/view/css/genome.css"],

// browsers Object. HASH ARRAY OF OPENED BROWSERS
browsers : [],

constructor : function(args) {	

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
	this.setViewProjectCombo();

	// SET COMBO LISTENERS
	setTimeout(
		function (thisObj) { thisObj.setFeatureProjectCombo(); },
		10,
		this
	);

	// SET VIEW COMBO ONKEYCHANGE LISTENER
	setTimeout(
		function (thisObj) { thisObj.setOnkeypressListener(); },
		1000,
		this
	);

	// LOAD BROWSER
	setTimeout(
		function (thisObj) { thisObj.loadBrowser(thisObj.getProject(), thisObj.getView()); },
		2000,
		this
	);

	var thisObject = this;
	this.mainTab.onClose = function() {

		var returned = confirm("Close this tab?");

		return returned;
	};
},

loadEval : function (filepath) {

	// ADD STAGE TO stage TABLE IN REMOTE DATABASE
	var url = Agua.htmlUrl + filepath;
	var timeout = 3000;

	// SEND TO SERVER
	dojo.xhrGet(
		{
			url: url,
			sync: true,
			handleAs: "text",
			load: function(response) {
				eval(response);
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);	
},


// GETTERS
getProject : function () {
	return this.viewProjectCombo.get('value');
},
getWorkflow : function () {
	return this.workflowCombo.get('value');
},
getView : function () {

	return this.viewCombo.get('value');
},
getViewFeature : function () {
	return this.featureList.get('value') ?
		this.featureList.get('value') : '' ;
},
getBuild : function () {

	return this.buildLabel.innerHTML;
},
getSpecies : function () {
	return this.speciesLabel.innerHTML;
},

getFeatureBuild : function () {

	var speciesBuild = this.speciesCombo.get('value');

	if ( speciesBuild.match(/^(\S+)\(([^\)]+)\)$/) )
		return  speciesBuild.match(/^(\S+)\(([^\)]+)\)$/)[2];
},

getFeatureSpecies : function () {

	var speciesBuild = this.speciesCombo.get('value');

	if ( speciesBuild.match(/^(\S+)\(([^\)]+)\)$/) )
		return speciesBuild.match(/^(\S+)\(([^\)]+)\)$/)[1];
},
getFeature : function () {
	return this.featureCombo.get('value');
},

getFeatureProject : function () {
	return this.featureProjectCombo.get('value');
},
getFeatureWorkflow : function () {
	return this.workflowCombo.get('value');
},

// VIEW METHODS
setOnkeypressListener : function () {

	// SET ONKEYPRESS LISTENER
	var thisObject = this;
	this.viewCombo._onKeyPress = function(event){

		// summary: handles keyboard events
		var key = event.charOrCode;			
		if ( key == 13 )
		{
			thisObject.workflowCombo._hideResultList();

			var projectName = thisObject.viewProjectCombo.get('value');
			var viewName = thisObject.viewCombo.get('value');

			// STOP PROPAGATION
			//event.stopPropagation();

			var isView = Agua.isView(projectName, viewName);

			if ( isView == false )	thisObject.addView(projectName, viewName);

			if ( thisObject.viewCombo._popupWidget != null )
			{
				thisObject.viewCombo._showResultList();
			}
		}

		// STOP PROPAGATION
		//event.stopPropagation();
	};
},

addView : function (projectName, viewName) {	
	var viewObject = new Object;
	viewObject.project = projectName;
	viewObject.view = viewName;
	var success = Agua.addView(viewObject);
	if ( success == true ) {
	}
	else {
		return;
	}

	// SET SPECIES LABEL TO BLANK
	this.speciesLabel.innerHTML = '';
	this.buildLabel.innerHTML = '';

	// ADD STAGE TO stage TABLE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "view.cgi";
	viewObject.username = Agua.cookie('username');
	viewObject.sessionId = Agua.cookie('sessionId');
	viewObject.mode = "addView";

	var callback = dojo.hitch(this, function() { this.setViewCombo(projectName, viewName)});

	this.doPut({ url: url, query: viewObject, callback: callback });	
},

removeView : function () {

	var projectName = this.getProject();
	var viewName = this.getView();

	// ADD TO DATA MODEL (LOCAL AND REMOTE)
	var viewObject = new Object;
	viewObject.project = projectName;
	viewObject.view = viewName;

	var success = Agua.removeView(viewObject);
	if ( success == true ) {
	}
	else {
		return;
	}

	// ADD STAGE TO stage TABLE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "view.cgi";
	viewObject.username = Agua.cookie('username');
	viewObject.sessionId = Agua.cookie('sessionId');
	viewObject.mode = "removeView";

	var callback = dojo.hitch(this, function() { this.setViewCombo(projectName)});

	this.doPut({ url: url, query: viewObject, callback: callback });	
},

setViewProjectCombo : function (projectName) {

	var projectNames = Agua.getProjectNames();
	if ( ! projectNames )
	{
		return;
	}

	// DO DATA ARRAY
	var data = {identifier: "name", items: []};
	for ( var i in projectNames )
		data.items[i] = { name: projectNames[i]	};
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );
	this.viewProjectCombo.store = store;	

	// SET PROJECT IF NOT DEFINED TO FIRST ENTRY IN projects
	if ( projectName == null || ! projectName)	projectName = projectNames[0];	
	this.viewProjectCombo.setValue(projectName);			

	// SET CSS
	this.viewProjectCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.viewProjectCombo.wrapperClass = "view dijitPopup";
	this.viewProjectCombo.itemHeight = 30;		

	if ( projectName == null )	projectName = this.viewProjectCombo.get('value');

	// RESET THE WORKFLOW COMBO
	this.setViewCombo(projectName);
},


setViewCombo : function (projectName, viewName) {

	this.killPopup(this.viewCombo);

	// DO COMBO WIDGET SETUP	
	this.inherited(arguments);

	// SET CSS
	this.viewCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.viewCombo.wrapperClass = "view dijitPopup";
	this.viewCombo.itemHeight = 30;		

	// SET VIEW NAME IF NOT DEFINED
	if ( viewName == null )
	{
		var views = Agua.getViewsByProject(projectName);
		if ( views.length > 0 ) viewName = views[0].view;
	}

	this.setSpeciesLabel(projectName, viewName);
},

setSpeciesLabel : function (projectName, viewName) {
	// SET SPECIES AND BUILD LABELS


	this.speciesLabel.innerHTML = '';
	this.buildLabel.innerHTML = '';

	var speciesHash = Agua.getViewSpecies(projectName, viewName);
	if ( speciesHash == null )
	{
		this.setFeatureList([]);
		return;
	}

	this.speciesLabel.innerHTML = speciesHash.species || '';
	this.buildLabel.innerHTML = speciesHash.build || '';

	// SET SPECIES COMBO VALUE
	if ( speciesHash.species != null && speciesHash.build )
	{
		var setValue = speciesHash.species + "(" + speciesHash.build + ")";
		this.speciesCombo.set('value', setValue);
	}

	// SET FEATURE LIST
	var viewfeatures = Agua.getViewFeatures(projectName, viewName);
	var featureNames = new Array;
	for ( var i = 0; i < viewfeatures.length; i++ )
		featureNames.push(viewfeatures[i].feature);

	this.setFeatureList(featureNames);
},

// FEATURE METHODS
setFeatureProjectCombo : function (projectName, workflowName) {

	var projectNames = Agua.getViewProjects();

	// DO DATA ARRAY
	var data = {identifier: "name", items: []};
	for ( var i in projectNames )
		data.items[i] = { name: projectNames[i]	};
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );
	this.featureProjectCombo.store = store;	

	// SET PROJECT IF NOT DEFINED TO FIRST ENTRY IN projects
	if ( projectName == null || ! projectName)	projectName = projectNames[0];	
	this.featureProjectCombo.setValue(projectName);			

	// SET CSS
	this.featureProjectCombo.popupClass = "feature featureCombo dijitReset dijitMenu";
	this.featureProjectCombo.wrapperClass = "feature dijitPopup";
	this.featureProjectCombo.itemHeight = 30;		

	if ( projectName == null )	projectName = this.featureProjectCombo.get('value');

	// RESET THE WORKFLOW COMBO
	this.setWorkflowCombo(projectName);
},

setWorkflowCombo : function (projectName, workflowName) {

	if ( projectName == null || ! projectName )
	{
		return;
	}

	// GET WORKFLOW COMBO
	if ( this.workflowCombo == null )
	{
		return;
	}

	// CREATE THE DATA FOR A STORE		
	var workflowNames = Agua.getViewProjectWorkflows(projectName);

	// RETURN IF workflowNames NOT DEFINED
	if ( ! workflowNames )
	{
		return;
	}		

	// CREATE store
	var data = {identifier: "name", items: []};
	for ( var i in workflowNames )
		data.items[i] = { name: workflowNames[i]	};
	var store = new dojo.data.ItemFileReadStore( { data: data } );
	this.workflowCombo.store = store;

	// START UP COMBO AND SET SELECTED VALUE TO FIRST ENTRY IN workflowNames IF NOT DEFINED 
	if ( workflowName == null || ! workflowName )	workflowName = workflowNames[0];

	this.workflowCombo.startup();
	this.workflowCombo.set('value', workflowName);			

	// SET CSS
	this.workflowCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.workflowCombo.wrapperClass = "view dijitPopup";
	this.workflowCombo.itemHeight = 30;		

	if ( projectName == null ) projectName = this.viewProjectCombo.get('value');
	if ( workflowName == null ) workflowName = this.workflowCombo.get('value');

	// RESET THE VIEW COMBO
	this.setSpeciesCombo(projectName, workflowName);
},

setSpeciesCombo : function (projectName, workflowName, speciesName, buildName) {

	// SET CSS
	this.viewCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.viewCombo.wrapperClass = "view dijitPopup";
	this.viewCombo.itemHeight = 30;		

	// SET DROP TARGET (LOAD MIDDLE PANE, BOTTOM)
	if ( projectName == null ) projectName = this.featureProjectCombo.get('value');
	if ( workflowName == null ) workflowName = this.workflowCombo.get('value');

	var viewfeatures = Agua.getViewWorkflowFeatures(projectName, workflowName);
	if ( viewfeatures == null || viewfeatures.length == 0 )
	{
		return;
	}

	// GET SPECIES+BUILD NAMES
	var speciesBuildNames = new Array;
	for ( var i = 0; i < viewfeatures.length; i++ )
	{
		speciesBuildNames.push(viewfeatures[i].species + "(" + viewfeatures[i].build + ")");
	}
	speciesBuildNames = this.uniqueValues(speciesBuildNames);

	// SET SPECIES+ BUILD NAME
	var speciesBuildName;
	if ( speciesName == null || ! speciesName
		|| buildName == null || ! buildName )
	{
		speciesBuildName = speciesBuildNames[0];
		speciesName = viewfeatures[0].species;
		buildName = viewfeatures[0].build;
	}
	else {
		speciesBuildName = speciesName + "(" + buildName + ")";
	}

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in speciesBuildNames )
	{
		data.items[i] = { name: speciesBuildNames[i]	};
	}
	var store = new dojo.data.ItemFileReadStore( { data: data } );
	this.speciesCombo.store = store;

	// START UP COMBO (?? NEEDED ??)
	this.speciesCombo.startup();
	this.speciesCombo.set('value', speciesBuildName);			

	this.setFeatureCombo(projectName, workflowName, speciesName, buildName);
},

setFeatureCombo : function (projectName, workflowName, speciesName, buildName) {

	// SET CSS
	this.featureCombo.popupClass = "view viewCombo dijitReset dijitMenu";
	this.featureCombo.wrapperClass = "view dijitPopup";
	this.featureCombo.itemHeight = 30;		

	if ( projectName == null || ! projectName 
		|| workflowName == null || ! workflowName 
		|| speciesName == null || ! speciesName 
		|| buildName == null || ! buildName )
	{
		return;
	}

	// SAVE THE SELECTED VALUE IF ONE EXISTS
	var featureName = this.featureCombo.get('value');

	// CREATE THE DATA FOR A STORE		
	var featureNames = 	Agua.getViewSpeciesFeatureNames(projectName, workflowName, speciesName, buildName);
	if ( ! featureNames )
	{
		return;
	}

	// CREATE store
	var data = {identifier: "name", items: []};
	for ( var i in featureNames )
		data.items[i] = { name: featureNames[i]	};
	var store = new dojo.data.ItemFileReadStore( { data: data } );
	this.featureCombo.store = store;

	// START UP COMBO AND SET SELECTED VALUE TO FIRST ENTRY IN featureNames IF NOT DEFINED 
	if ( featureName == null || ! featureName )	featureName = featureNames[0];

	this.featureCombo.startup();
	this.featureCombo.set('value', featureName);			
},


setFeatureList : function (featureNames) {

	// SET CSS
	this.featureList.popupClass = "view viewCombo dijitReset dijitMenu";
	this.featureList.wrapperClass = "view dijitPopup";
	this.featureList.itemHeight = 30;		

	var data = {identifier: "name", items: []};
	if ( featureNames == null || featureNames.length == 0 )
	{
	}
	else{
		for ( var i in featureNames )
			data.items[i] = { name: featureNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( { data: data } );
	this.featureList.store = store;

	// START UP COMBO AND SET SELECTED VALUE TO FIRST ENTRY 
	this.featureList.startup();
	this.featureList.set('value', featureNames[0]);			

},

removeViewFeature : function () {

	// ADD TO DATA MODEL (LOCAL AND REMOTE)
	var featureObject = new Object;
	featureObject.project = this.getProject();
	featureObject.view = this.getView();
	featureObject.feature = this.getViewFeature();
	featureObject.species = this.getFeatureSpecies();
	featureObject.build = this.getFeatureBuild();

	Agua.removeViewFeature(featureObject);

	this.setViewCombo(this.getProject());
},

addViewFeature : function () {

	// ADD TO LOCAL DATA MODEL ON CLIENT
	var featureName = this.featureCombo.get('value');

	var featureObject = new Object;
	featureObject.sourceproject = this.getFeatureProject();
	featureObject.sourceworkflow = this.getFeatureWorkflow();
	featureObject.species = this.getFeatureSpecies();
	featureObject.build = this.getFeatureBuild();
	featureObject.feature = this.getFeature();
	featureObject.project = this.getProject();
	featureObject.view = this.getView();

	// ADD the FEATURE TO THE LOCAL DATA MODEL ON THE CLIENT
	var success = Agua.addViewFeature(featureObject);
	if ( success == true ) {
	}
	else {
		return;
	}

	// ADD TO viewfeature TABLE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "view.cgi";
	featureObject.username = Agua.cookie('username');
	featureObject.sessionId = Agua.cookie('sessionId');
	featureObject.mode = "addViewFeature";

	var callback = dojo.hitch(this, function() { this.setViewCombo(featureObject.project, featureObject.view)});

	this.doPut({ url: url, query: featureObject, callback: callback });	
},


// BROWSER METHODS
loadBrowser : function (projectName, viewName) {

	if ( projectName == null )	projectName = this.getProject();
	if ( viewName == null )		viewName = this.getView();

	// IF SPECIES OR BUILD ARE NOT DEFINED THEN THIS VIEW
	// DOES NOT HAVE ANY FEATURES ADDED TO IT YET
	if ( ! this.getSpecies() )
	{
		return;
	}
	if ( ! this.getBuild() )
	{
		return;
	}

	var username = Agua.cookie('username');
	var refseqfile = "plugins/view/jbrowse/users"
					+ "/" + username
					+ "/" + projectName
					+ "/" + viewName
					+ "/data/refSeqs.js";
	var trackinfofile = "plugins/view/jbrowse/users"
					+ "/" + username
					+ "/" + projectName
					+ "/" + viewName
					+ "/data/trackInfo.js";

	// LOAD refSeqs AND trackInfo JSON FILES
	this.loadEval(trackinfofile);
	this.loadEval(refseqfile);
	//this.loadEval("plugins/view/jbrowse/species/human/hg19/data/trackInfo.js");
	//this.loadEval("plugins/view/jbrowse/species/human/hg19/data/refSeqs.js");

	// CHECK INPUTS
	if ( projectName == null || viewName == null )
	{
		return;
	}

	var viewObject = Agua.getViewObject(projectName, viewName);
	var location	=	viewObject.chromosome + ":" + viewObject.start + "..." + viewObject.stop;
	var trackList	=	viewObject.tracklist;

	// GET UNIQUE ID FOR THIS MENU TO BE USED IN DND SOURCE LATER
	var objectName = "plugins.view.View.jbrowse.Browser";
	var browserId = dijit.getUniqueId(objectName.replace(/\./g,"_"));

	var dataRoot = "plugins/viewName/jbrowse/" + viewObject.species + "/";

	var speciesName = this.getSpecies();
	var buildName = this.getBuild();
	var username = Agua.cookie('username');

	var b = new plugins.view.jbrowse.js.Browser({
		parentWidget : this,
		viewObject : viewObject,
		speciesName: speciesName,
		buildName  : buildName,
		species    : speciesName,
		build      : buildName,
		refSeqs: refSeqs,
		trackData: trackInfo,
		baseUrl : "plugins/view/jbrowse/",
		//dataRoot : "plugins/view/jbrowse/species/" + viewObject.species + "/" + viewObject.build + "/",
		dataRoot : "plugins/view/jbrowse/users/" + username + "/" + viewObject.project + "/" + viewObject.view + "/",
		browserRoot : "plugins/view/jbrowse/",
		//defaultLocation : "chr2:10000000..100000000",
		//trackList : "Affy_Exon_Probes,Microsatellite,Encode_Common_CNV",
		//trackList : "Chromosome_Band,OMIM_Genes,UCSC_Genes",
		//trackList : "vegaGene,tRNAs,gap,ccdsGene",
		trackList : trackList,
		//showTracks : "Assembly,Gap",

		attachWidget : this.rightPane
		//dataRoot : "plugins/viewName/jbrowse/users/syoung/Project1/Workflow1/"
	});

	// ADD TO this.browsers ARRAY		
	//this.addBrowser(b, projectName, viewName);

	b.showTracks("Assembly","Gap");

	b.navigateTo("chr2:10000000-20000000");

	if ( location != null && location != '' )
	{
		b.navigateTo(location);
	}

}, // 	loadBrowser 



deleteView : function () {

	var project = this.getProject();
	var workflow = this.getWorkflow();
	var view = this.getView();



},

// FIRE COMBO HANDLERS
fireViewProjectCombo : function() {
	var projectName = this.viewProjectCombo.get('value');
	this.setViewCombo(projectName);
},

fireViewCombo : function () {
// ONCHANGE IN VIEW COMBO FIRED
	if ( ! this.viewComboFired == true )
	{
		this.viewComboFired = true;
	}
	else {

		var projectName = this.getProject();
		var viewName = this.getView();
		this.setSpeciesLabel(projectName, viewName);

		this.loadBrowser(projectName, viewName);
	}
},

fireFeatureProjectCombo : function() {
	if ( ! this.featureProjectComboFired == true )
	{
		this.featureProjectComboFired = true;
	}
	else {
		var projectName = this.featureProjectCombo.get('value');
		this.setWorkflowCombo(projectName);
	}
},

fireWorkflowCombo : function() {
	if ( ! this.workflowComboFired == true )
	{
		this.workflowComboFired = true;
	}
	else {
		var projectName = this.featureProjectCombo.get('value');
		var workflowName = this.workflowCombo.get('value');
		this.setSpeciesCombo(projectName, workflowName);
	}
},

fireSpeciesCombo : function () {
	if ( ! this.speciesComboFired == true )
	{
		this.speciesComboFired = true;
	}
	else {
		var projectName = this.viewProjectCombo.get('value');
		var workflowName = this.workflowCombo.get('value');
		var speciesName = this.getSpecies();
		var buildName = this.getBuild();

		if ( speciesName == null || buildName == null )
		{
			return;
		}

		this.setFeatureCombo(projectName, workflowName, speciesName, buildName);
	}
}




}); // end of plugins.view.View
