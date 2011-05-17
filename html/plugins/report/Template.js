
dojo.provide("plugins.report.Template");

// REFERENCE SITES

// Object Oriented concepts and inheritance
// http://dojotoolkit.org/book/dojo-book-0-4/part-3-dojo-programming-model/object-oriented-concepts-and-inheritance

// BIOMART
//http://www.ensembl.org/biomart/martview/e61a39e4e8e306354f8f2a7a70bbc53c/e61a39e4e8e306354f8f2a7a70bbc53c

// DOJO DEMOS
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/demos/
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/dojox/form/tests/test_SelectStack.html

//Real World Dojo part 5: Custom Components
//http://www.dotnetmafia.com/blogs/jamesashley/archive/2008/10/28/761.aspx

dojo.declare(
    "plugins.report.Template",
    [dijit._Widget, dijit._Templated],
    {
        //Path to the template of this widget. dijit._Templated uses this to
        //snatch the template from the named file via a synchronous call.
        // This provides the relative path to the template file for the widget. Note that
        // fetching the template for a widget requires a synchronous network call, although
        // Dojo will cache the template string after it is fetched (unless a custom build
        // with interned template strings is used).
        templatePath: null,

        // Calls dijit._Templated.widgetsInTemplate
        // If widgets are defined inside of the template (path or string), this
        // value should be specified and set to true so that the Dojo parser will
        // find and instantiate those widgets. This value is false by default.
        // While this property is easy to overlook, know that defining widgets
        // inside of other widget?s templates can often be very useful.
        widgetsInTemplate : true,

        // SLOTS TO BE POPULATED BY THE DESCENDENT OF THIS CLASS
        title : '',
        id : '',
		filename: '',


		// STORE ELEMENTS OF INTEREST
		elements : [],
		elementHash: null,
        elementObjects : [],
		elementValues : null,

		// WIDGET TYPES FOR WHICH THE DOJO CLASSES ARE KNOWN
		knownTypes : null,

		// HASH OF FUNCTIONS TO GET AND SET VALUES IN WIDGETS AND DOM NODES
		valueFunctions : null,

		// CSS FILES
        cssFiles : [],

        //Any initialization code would go here in the constructor. dijit._Widget and
        //dijit._Templated do not have parameters in their constructors, so
        //there wouldn't be any multiple-inheritance complications
        //if you were to include some paramters here.
        constructor : function (args)
        {

            this.id         = args.id;
			this.filename   = args.filename;
			//this.url        = args.url;

            // LOAD CSS FOR WIDGETS
            this.loadCSS();
        },


        //Inherited from dijit._Widget and called just before template
        //instantiation in buildRendering. This method is especially useful
        //for manipulating the template before it becomes visible.
        postMixInProperties: function()
        {
        },


        //You can override this method to manipulate widget once it is
        //placed in the UI, but be warned that any child widgets contained
        //in it may not be ready yet.        
        postCreate: function()
		{
		},


        //Called after the widget's children and all other widgets on the
        //page have been created. Provides an opportunity to manipulate child
        //widgets before they're displayed.
        // NOW THAT ALL THE CHILD WIDGETS HAVE BEEN CREATED, ADD EACH CHILD WIDGET
        // TO ITS RELATED widgetObject IN THE GLOBAL this.elementObjects HASH

		// USE this.inherited TO CALL THIS AT THE START OF THE SUBCLASS OBJECT'S
		// startup METHOD THEN SET UP THE PAGE WITH MORE COMPLICATED LOGIC 

        startup : function ()
        {

			// SET ELEMENT HASH
			this.setElementHash();

            // BUILD HASH OF ELEMENTS (name: ELEMENT OBJECT) AND SET
			// UNIQUE IDS FOR EACH 
            this.setElementObjects();

			// SET ELEMENTS: ADD SPECIFIC PROPERTIES, BASIC LOGIC AND
			// VALUES TO WIDGETS OR NODES
			this.setElements();
        },


		// setElements
		//
		// ADD SPECIFIC PROPERTIES AND VALUES TO WIDGETS OR NODES
		// 
		// USES THE FACT THAT this[functionName] CAN BE A FUNCTION

		setElements : function ()
		{

			for ( name in this.elementHash )
			{

				if ( this[name] )
				{
					// RUN THE FUNCTION
					this[name](this.elementObjects[name], name);

					// **** DEBUG : SHOW THE VALUE 
					if ( this.elementObjects[name].valueFunction )
					{
					}
				}
				else
				{
				}
			}
		},


		// SET HASH OF ELEMENTS OF INTEREST WITH IDS IN TEMPLATE.HTML
		setElementHash : function ()
		{

			// SET HASH OF ELEMENT IDS
			this.elementHash = new Object;
			for ( var i = 0; i < this.elements.length; i++ )
			{
				if ( this.elementHash[this.elements[i]] )
				{
				}
				else
				{
					this.elementHash[this.elements[i]] = 1;
				}
			}
		},


		// setElementObjects
		//
		// BUILD HASH OF ELEMENTS (name: ELEMENT OBJECT) AND SET
		//
		// UNIQUE IDS FOR EACH 
		//
		// NOTES:
		// dojo.query AND forEach
		// http://www.dojotoolkit.org/book/dojo-book-0-9/part-3-programmatic-dijit-and-dojo/functions-used-everywhere/dojo-foreach

        setElementObjects : function ()
        {
			var thisTemplate = this;

			dojo.query("*", this.containerNode).forEach(

				function(element) {

					if ( element.id && thisTemplate.elementHash[element.id] ) 
					{
						var widget = dijit.byId(element.id);
						if ( widget != null )
						{
							thisTemplate.elementObjects[element.id] = widget;
						}
						else
						{
							thisTemplate.elementObjects[element.id] = element;
						}


						// SET UNIQUE ID
						var uniqueId = dijit.getUniqueId(thisTemplate.id);
						element.id = uniqueId;
					}
				}
			);

			return;
        },




		// loadCSS
		//
		// LOAD CSSFILES IN this.cssFiles ARRAY
		//

        loadCSS : function()
        {

            // LOAD CSS
            var cssFiles = this.cssFiles;

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
        },




		// loadWidgetValues
		//
		// LOAD A LIST OF name : value HASHED VALUES INTO ELEMENTS IN THE TEMPLATE
        //

        loadWidgetValues : function (response, ioArgs)
        {

			for ( var name in response )
			{

                var value = response[name];
				if ( this.elementObjects[name] )
				{
					this.elementObjects[name].value(this.elementObjects[name], value);
				}
				else
				{
				}
				//this.elementObjects[name].value(this.elementObjects[name], "'" + value + "'");
			}
        },


		valueFunction : function (name)
		{
			// SET VALUE FUNCTION FOR saveReport LATER
			switch (name)
			{
				case "combobox" : case "editor" : case "checkbox" :

					return function(value) {
						if ( value )
						{
							return this.setValue(value);
						}
						else
						{
							var returnValue = this.getValue();
							if ( returnValue != null && returnValue != '' )
							{
								returnValue = returnValue.replace(/"/g, "'");
							}
							return returnValue;
						}
					};

				case "null" :
					return function(value) { return null; };

				case "radio" : return function(value) {
						return value ? this.checked = value : this.checked;
					};

				case "textinput" : return function(value) {
						return value ? this.value = value : this.value;
					};

				case "div" : return function(value) {
						return value ? this.innerHTML = value : this.innerHTML;
					};			
			}

			return 0;
		}
    }


); // plugins.report.Template





//		// method over-ride
//		buildRendering: function(){
//			// summary:
//			//		Construct the UI for this widget from a template, setting this.domNode.
//
//			// Lookup cached version of template, and download to cache if it
//			// isn't there already.  Returns either a DomNode or a string, depending on
//			// whether or not the template contains ${foo} replacement parameters.
//
//
//			var cached = dijit._Templated.getCachedTemplate(this.templatePath, this.templateString, this._skipNodeCache);
//
//
//			var node;
//			if(dojo.isString(cached)){
//
//
//				node = dijit._Templated._createNodesFromText(this._stringRepl(cached))[0];
//			}else{
//				// if it's a node, all we have to do is clone it
//				node = cached.cloneNode(true);
//			}
//
//
//			this.domNode = node;
//
//			// recurse through the node, looking for, and attaching to, our
//			// attachment points and events, which should be defined on the template node.
//			this._attachTemplateNodes(node);
//
//
//
//
//			var source = this.srcNodeRef;
//
//
//			if(source && source.parentNode){
//
//
//				source.parentNode.replaceChild(node, source);
//			}
//
//
//
//// IF IDS HAVE ALREADY BY SET BY setElementObjects, THIS GENERATES THIS ERROR:
//// Tried to register widget with id==snpReport1_27 but that id is already registered
//
//			if(this.widgetsInTemplate){
//
//
//				var cw = (this._supportingWidgets = dojo.parser.parse(node));
//				this._attachTemplateNodes(cw, function(n,p){
//					return n[p];
//				});
//			}
//
//
//			this._fillContent(source);
//		},
//
//
//
//
//		_attachTemplateNodes: function(rootNode, getAttrFunc){
//			// summary: Iterate through the template and attach functions and nodes accordingly.	
//			// description:		
//			//		Map widget properties and functions to the handlers specified in
//			//		the dom node and it's descendants. This function iterates over all
//			//		nodes and looks for these properties:
//			//			* dojoAttachPoint
//			//			* dojoAttachEvent	
//			//			* waiRole
//			//			* waiState
//			// rootNode: DomNode|Array[Widgets]
//			//		the node to search for properties. All children will be searched.
//			// getAttrFunc: function?
//			//		a function which will be used to obtain property for a given
//			//		DomNode/Widget
//
//
////			console.dir(this);
//
//
//			getAttrFunc = getAttrFunc || function(n,p){ return n.getAttribute(p); };
//
//			var nodes = dojo.isArray(rootNode) ? rootNode : (rootNode.all || rootNode.getElementsByTagName("*"));
//			var x = dojo.isArray(rootNode) ? 0 : -1;
//			var attrs = {};
//			for(; x<nodes.length; x++){
//				var baseNode = (x == -1) ? rootNode : nodes[x];
//				if(this.widgetsInTemplate && getAttrFunc(baseNode, "dojoType")){
//					continue;
//				}
//				// Process dojoAttachPoint
//				var attachPoint = getAttrFunc(baseNode, "dojoAttachPoint");
//				if(attachPoint){
//					var point, points = attachPoint.split(/\s*,\s*/);
//					while((point = points.shift())){
//						if(dojo.isArray(this[point])){
//							this[point].push(baseNode);
//						}else{
//							this[point]=baseNode;
//						}
//					}
//				}
//
//
//				// Process dojoAttachEvent
//				var attachEvent = getAttrFunc(baseNode, "dojoAttachEvent");
//				if(attachEvent){
//					// NOTE: we want to support attributes that have the form
//					// "domEvent: nativeEvent; ..."
//					var event, events = attachEvent.split(/\s*,\s*/);
//					var trim = dojo.trim;
//					while((event = events.shift())){
//						if(event){
//							var thisFunc = null;
//							if(event.indexOf(":") != -1){
//								// oh, if only JS had tuple assignment
//								var funcNameArr = event.split(":");
//								event = trim(funcNameArr[0]);
//								thisFunc = trim(funcNameArr[1]);
//							}else{
//								event = trim(event);
//							}
//							if(!thisFunc){
//								thisFunc = event;
//							}
//							this.connect(baseNode, event, thisFunc);
//						}
//					}
//				}
//
//
//				// waiRole, waiState
//				var role = getAttrFunc(baseNode, "waiRole");
//				if(role){
//					dijit.setWaiRole(baseNode, role);
//				}
//				var values = getAttrFunc(baseNode, "waiState");
//				if(values){
//					dojo.forEach(values.split(/\s*,\s*/), function(stateValue){
//						if(stateValue.indexOf('-') != -1){
//							var pair = stateValue.split('-');
//							dijit.setWaiState(baseNode, pair[0], pair[1]);
//						}
//					});
//				}
//			}
//
//
//		},
//


		// DON'T DO THIS HERE --> LET THE USER SET THE VALUE FUNCTION FOR EACH ELEMENT
		// IN THE this[elementName] FUNCTION

		////////setKnownTypes : function()
		////////{
		////////	this.knownTypes = new Array;
		////////	knownTypes["checkbox"] = ["dijit.form.CheckBox"];
		////////
		////////},
		////////
		////////// GET ELEMENT NAME:
		//////////
		//////////	INPUT 	:	WIDGET TYPE (E.G., checkbox)
		//////////
		////////// 	OUTPUT	:	WIDGET CLASS NAME (E.G., dijit.form.CheckBox)
		//////////
		//////////
		////////getElementName : function (type)
		////////{
		////////
		////////	
		////////
		////////},




        // EXAMPLE OF xhrPost
        // http://www.dojoforum.com/2007/10/11/dojo-example-xhrget-and-xhrpost

//////		// STUB FOR LOADING REPORT DATA FROM A WORKFLOW
//////		loadReport : function()
//////		{
////////            	
////////			// GENERATE QUERY
////////			var query = new Object;
////////			query.username = Agua.cookie('username');
////////			query.sessionId = Agua.cookie('sessionId');
////////			query.project = Agua.cookie('project');
////////			query.workflow = Agua.cookie('workflow');
////////			query.report = Agua.cookie('Report1');
////////			query.mode = "loadReport";
////////			query.module = "Report::SNP";
////////	
////////			// POST DATA TO SERVER
////////            var template = this;
////////			dojo.xhrPost(
////////				{
////////					//url: this.url,
////////                    url : "plugins/report/templates/SNP.json",
////////					contentType: "text",
////////					handleAs: "json",
////////					postData: dojo.toJson(query),
////////					//handleAs: "json-comment-filtered",
////////					timeout: 2000,
////////					load: function(response, ioArgs) {
////////                        template.loadWidgetValues(response, ioArgs);
////////                    },
////////					error: function(response, ioArgs) {
////////						return response;
////////					}
////////				}
////////			);
//////			
//////		},


