
dojo.provide("plugins.report.Grid");

dojo.require("dijit.dijit"); // optimize: load dijit layer

// TITLE PANE
dojo.require("dijit.TitlePane");

// GRID
dojo.require("dojo.data.ItemFileWriteStore");
dojo.require("dojox.grid.DataGrid");

dojo.require("dijit.Tooltip");
dojo.require("dijit.Menu");

// DEPTH NUMBER SPINNER
dojo.require("dijit.form.NumberSpinner");

// PARSER
dojo.require("dojo.parser");

// INHERITED CLASSES
dojo.require("plugins.core.Common");

dojo.declare(
    "plugins.report.Grid",
    [ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
    //Path to the template of this widget. 
    templatePath: dojo.moduleUrl("plugins", "report/templates/grid.html"),

    // Calls dijit._Templated.widgetsInTemplate
    widgetsInTemplate : true,

    // OR USE @import IN HTML TEMPLATE
    cssFiles : [ "dojo-1.5.0/dojox/grid/resources/Grid.css", "dojo-1.5.0/dojox/grid/resources/tundraGrid.css", "plugins/report/css/grid.css" ],

	// MAX. NO. OF USERS TO DISPLAY AT ANY ONE TIME IN DRAG SOURCE
	maxResults : 40,

    // Any initialization code would go here in the constructor.
    // plugins.report.Template and its superclasses dijit._Widget and
    // dijit._Templated do not have parameters in their constructors, so
    // there wouldn't be any multiple-inheritance complications
    // if you were to include some paramters here.
    constructor : function (args)
    {
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

        this.startup();
    },

    // startup
    //
    //
    startup : function ()
    {

        // SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
        this.inherited(arguments);

        // LOAD CSS STYLE FILES
        this.loadCSS();

        // ADD TO TAB CONTAINER		
        this.attachNode.appendChild(this.containerNode);

        this.setResultsSlider();
        this.setResultsSpinner();
        this.showMaxResults();
    },

    // INITIALISE SLIDER TO SELECT BOUNDARIES OF RESULTS RANGE
    setResultsSlider : function ()
    {

        // ONMOUSEUP
        var thisObject = this;
        dojo.connect(this.gridSlider, "onMouseUp", dojo.hitch(this, function(e)
        {
            var position = parseInt(this.gridSlider.getValue());
            thisObject.displayResults(position);
        }));
    },

	// INITIALISE SPINNER TO SELECT WIDTH OF RESULTS RANGE
	setResultsSpinner : function ()
	{
		this.gridSpinnerNode.id = '';
		this.gridSpinner = new dijit.form.NumberSpinner(
			{
                value: this.maxResults,
                smallDelta: "10",
				constraints: { min: 20, max: 500, places:0},
				size : 1
			},
			this.gridSpinnerNode
		);
		this.gridSpinner.domNode.setAttribute('style', "width: 50px;");

		//this.gridSpinner.downArrowNode.setAttribute('style', "background: transparent url(plugins/report/images/spriteArrows.png) no-repeat scroll !important; width: 15px;"); 

		//this.gridSpinner.downArrowNode.setAttribute('style', "background: transparent url(http://localhost/agua/0.4/plugins/view/images/spriteArrows.png) no-repeat scroll !important; width:17px;"); 
		//this.gridSpinner.downArrowNode.setAttribute('class', 'downArrow'); 

        dojo.connect(this.gridSpinner, "onBlur", dojo.hitch(this, function ()
        {
            this.maxResults = this.gridSpinner.value;

            var position = this.resultStart.innerHTML;

            this.displayResults(position);
        }));


	},

	// INITIALISE SPINNER TO SELECT WIDTH OF RESULTS RANGE
	updateResultsSpinner : function ()
	{
		if ( this.gridSpinner == null )
        {
            return;
        }

		//this.gridSpinner.domNode.setAttribute('style', "width: 50px;");
	},


    // SET INNERHTML OF this.showResults NODE TO this.maxResults
    showMaxResults : function ()
    {

        //this.showResults.innerHTML = this.maxResults;
    },


    // SET this.maxResults TO AN INTEGER
    setMaxResults : function (integer)
    {
        if ( integer == null || ! integer )
        {
            return;
        }

        integer = parseInt(integer);

        if ( integer < 1 )
        {
            return;
        }

        this.maxResults = integer;
    },

    // SET INNERHTML OF this.totalResults NODE TO AN INTEGER
    showTotalResults : function (integer)
    {

        this.totalResults.innerHTML = integer;
    },

    // SET INNERHTML OF this.resultStart NODE TO AN INTEGER
    showResultStart : function (integer)
    {

        this.resultStart.innerHTML = integer;
    },


    // SET INNERHTML OF this.resultStop NODE TO AN INTEGER
    showResultStop : function (integer)
    {

        this.resultStop.innerHTML = integer;
    },


    // LOAD DATA INTO this.data AND CALL displayResults
    load : function(data)
    {

        this.data = data;

        // DISPLAY RESULTS
        this.displayResults(null);
    },




	displayResults : function (position)
	{

        if ( this.doingDisplay == true )
        {
            return;
        }
        this.doingDisplay = true;

        // DISABLE GRID SLIDER
        this.gridSlider.attr('disabled', true);


		// REMOVE ALL EXISTING CONTENT
        // REMOVE EXISTING GRID NODE FROM GRID CONTAINER
        var gridContainer = this.outputResult;
        while(gridContainer.firstChild)
        {
           gridContainer.removeChild(gridContainer.firstChild);
        }

        // SET GRID ID AND REMOVE ANY EXISTING GRID
        var gridId = this.gridId;
        if ( dijit.byId(gridId) )
        {
            dijit.byId(gridId).destroy();
        }


		// SET position IF NOT DEFINED
		if ( position == null )	position = 0;


		// RETURN IF this.data IS NULL OR EMPTY
		if ( this.data == null || this.data.length == 0 )
		{
            this.doingDisplay = false;

            // ENABLE GRID SLIDER
            this.gridSlider.attr('disabled', false);

			return;
		}

		var MAXRESULTS = this.maxResults;
		var MULTIPLE = ( this.data.items.length - MAXRESULTS ) / 100;
		MULTIPLE = MULTIPLE ? MULTIPLE : 1;
		var start = parseInt(position * MULTIPLE);
        if ( start < 0 ) start = 0;
		var stop = parseInt( (position * MULTIPLE) + MAXRESULTS );
		if ( ! stop || stop > this.data.items.length )	stop = this.data.items.length;

        // DISPLAY START AND STOP POINTS
        this.showTotalResults(this.data.items.length);
        this.showResultStart(start);
        this.showResultStop(stop);


		// GENERATE USER DATA TO INSERT INTO DND USER TABLE
        var tempData = new Object;
        tempData.identifier = this.data.identifier;
        tempData.label = this.data.label;

		var itemsArray = new Array;
        for ( var j = start; j < stop; j++ )
		{
			itemsArray.push(this.data.items[j]);
		}
        tempData.items = itemsArray;



        // SEE EXAMPLE:    
        //http://localhost:8080/Bioptic0.2.5/html/tests/report-snp3.html

        // SET STORE USING DATA
        var gridStore = new dojo.data.ItemFileWriteStore({data: tempData});


        var layout = [
            {name: 'Entry', field: 'id', width: "30px"},
            {name: 'Name', field: 'name', width: "100px"},
            {name: 'Chr', field: 'chromosome', width: "40px"},
            {name: 'CCDS start', field: 'ccdsstart', width: "40px"},
            {name: 'CCDS stop', field: 'ccdsstop', width: "40px"},
            {name: 'Ref nt', field: 'referencenucleotide', width: "20px"},
            {name: 'Var nt', field: 'variantnucleotide', width: "20px"},
            {name: 'Depth', field: 'depth', width: "30px"},
            {name: 'Var freq', field: 'variantfrequency', width: "30px"},
            {name: 'Chr start', field: 'chromosomestart', width: "40px"},
            {name: 'Chr stop', field: 'chromosomestop', width: "40px"},
            {name: 'Sense', field: 'sense', width: "40px"},
            {name: 'Ref codon', field: 'referencecodon', width: "35px"},
            {name: 'Var codon', field: 'variantcodon', width: "35px"},
            {name: 'Ref aa', field: 'referenceaa', width: "50px"},
            {name: 'Var aa', field: 'variantaa', width: "50px"},
            {name: 'Strand', field: 'strand', width: "35px"},
            {name: 'SNP', field: 'snp', width: "40px"},
            {name: 'Score', field: 'score', width: "40px"},
            {name: 'dbSNP strand', field: 'dbsnpstrand', width: "40px"}
        ];

        var gridId = dojo.dnd.getUniqueId();

        grid = new dojox.grid.DataGrid({
            id: gridId,
            autoHeight: true,
            store: gridStore,
            structure: layout
        }, document.createElement('div'));

        // IMPORTANT: CLASS IS REQUIRED FOR PROPER DISPLAY
        dojo.addClass(grid.viewsHeaderNode, 'viewsHeaderNode');
        dojo.addClass(grid.domNode, 'reportSNPGrid');
        dojo.addClass(gridContainer, 'work-pane-split');


        grid.onCellContextMenu = function(e) {
                cellNode = e.cellNode;
        };
        grid.onHeaderContextMenu = function(e) {
            cellNode = e.cellNode;
        };    

//
//            var dojoxGridViewNode = grid.viewsNode.firstChild;
//			//dojoxGridViewNode children
//			//0
//			//	input.dojoxGridHiddenFocus on
//			//1
//			//	input.dojoxGridHiddenFocus on
//			//2
//			//	div.dojoxGridScrollbox
//            var dojoxGridScrollbox = dojoxGridViewNode.childNodes[2];


        //dojo.addClass(grid.headerContentNode, 'headerContentNode');



        // RESPONDS TO CLICKS ON THE GRID -- BUT CAN'T FIND GRID...
        this.grid = grid;

        var snpReport = this;
        dojo.connect(grid,'onClick', function getGridRow (n)
            {
                var grid = snpReport.grid;

                var row = grid.selection.getSelected()[0];

                var datarow = grid.getItem(n.rowIndex);


                var name = grid.store.getValue(datarow,"name");
                var chromosome = grid.store.getValue(datarow,"chromosome");
                var chromosomeStart = grid.store.getValue(datarow,"chromosomestart");
                var chromosomeStop = grid.store.getValue(datarow,"chromosomestop");

                //var t = grid.store.getValue(row,"Title");
                //console.debug("jjjj jjjj f: " +f+" "+t);

                return;
            }
        );

//            // ADD CONTEXT MENU
        var gridMenu = this.createGridMenu(grid);
        gridMenu.bindDomNode(grid.domNode);

        // prevent grid methods from killing the context menu event by implementing our own handler
        grid.onCellContextMenu = function(e) {

            cellNode = e.cellNode;
            console.dir(cellNode);
        };

        grid.onHeaderContextMenu = function(e) {

            cellNode = e.cellNode;
            console.dir(cellNode);
        };


        // APPEND GRID TO GRID CONTAINER AND STARTUP
        gridContainer.appendChild(grid.domNode);
        grid.startup();


        // ENABLE GRID SLIDER
        this.gridSlider.attr('disabled', false);

        this.doingDisplay = false;

    },   // load



    // ADD PROGRAMMATIC CONTEXT MENU
    createGridMenu : function (grid)
    {
        var menuId = dojo.dnd.getUniqueId();
        var dynamicMenu = new dijit.Menu( { id: menuId } );

        // ADD MENU TITLE
        dynamicMenu.addChild(new dijit.MenuItem( { label:"Options", disabled:false} ));
        dynamicMenu.addChild(new dijit.MenuSeparator());

        //// ONE OF FOUR WAYS TO DO MENU CALLBACK WITH ACCESS TO THE MENU ITEM AND THE CURRENT TARGET 	
        // 4. dojo.connect CALL
        //	REQUIRES:
        //		ADDED menu.currentTarget SLOT TO dijit.menu
        var mItem1 = new dijit.MenuItem(
            {
                label: "View",
                disabled: false
            }
        );
        dynamicMenu.addChild(mItem1);


        var snpReport =this;
        dojo.connect(mItem1, "onClick", function(e)
            {
                if( cellNode )
                {
                    var rowNode = cellNode.parentNode;
                    var tableNode = cellNode.parentNode.parentNode.parentNode;
                    var divNode = tableNode.parentNode;
                    var containerNode = tableNode.parentNode.parentNode;

                    var rowIndex;
                    for ( rowIndex = 0; rowIndex < containerNode.childNodes.length; rowIndex++ )
                    {
                        if ( containerNode.childNodes[rowIndex] == divNode )
                        {
                            break;
                        }
                    }

                    var datarow = grid.getItem(rowIndex);
                    var name = grid.store.getValue(datarow,"name");
                    var chromosome = grid.store.getValue(datarow,"chromosome");
                    var chromosomeStart = grid.store.getValue(datarow,"chromosomestart");
                    var chromosomeStop = grid.store.getValue(datarow,"chromosomestop");


                    snpReport.viewRow(rowIndex, grid);
                }
            }
        );

        // SEPARATOR
        dynamicMenu.addChild(new dijit.MenuSeparator());

        //	ADD run MENU ITEM
        var mItem2 = new dijit.MenuItem(
            {
                label: "Edit",
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

    },   // createGridMenu


    // DISPLAY THE GRID ITEM IN THE 'View' TAB
    viewRow : function (rowIndex, grid)
    {
        //var value = row[0];

        //var row = grid.selection.getSelected()[0];
        var datarow = grid.getItem(rowIndex);
        var name = grid.store.getValue(datarow,"name");
        var chromosome = grid.store.getValue(datarow,"chromosome");
        var chromosomeStart = grid.store.getValue(datarow,"chromosomestart");
        var chromosomeStop = grid.store.getValue(datarow,"chromosomestart");


        //var t = grid.store.getValue(row,"Title");
        //console.debug(f+" "+t);

        //var store = grid.store;

        //var row = store.fetch(
        //{
        //    query : { id : rowIndex }, 
        //    onComplete: function(items, request)
        //    {
        //        
        //        return 1;
        //    }
        //});

        //var row = grid.getItem(rowIndex);
        //
        //var value = row[0];
        //
        //var row = grid.selection.getSelected()[rowIndex];
        //console.debug("row:" + row);

        //var value = grid.model.data[row];
        //console.debug("value: " + value);

        //var erName = value.ER_Name;
        //console.debug("erName : " + erName);


    },

    // GRID REPORT GOES HERE
    output : function (widget, name)
    {
        var snpReport = this;

        // SET VALUE FUNCTION
        if ( ! this.elementObjects[name] )
        {
        }



        this.elementObjects[name].valueFunction = function(widgetObject, value)
        {
            if ( widgetObject.widget )
            {
                return value ? widgetObject.widget.setValue(value) : widgetObject.widget.getValue();
            }
            else if ( dojo.byId(widgetObject.id ) && ( value == null || ! value ) )
            {
                return dojo.byId(widgetObject.id).value; 
            }
        }
    }


}); // plugins.report.Grid



//setMenu : function ()
//{
//    var showTooltip = function(e)
//    {
//        if(gridTooltipEnabled)
//        {
//            var msg = "This is cell " + e.rowIndex + ", " + e.cellIndex;
//            dijit.showTooltip(msg, e.cellNode);
//        }
//    },
//
//    hideTooltip = function(e) {
//        dijit.hideTooltip(e.cellNode);
//        // FIXME: make sure that pesky tooltip doesn't reappear!
//        // would be nice if there were a way to hide tooltip without regard to aroundNode.
//        dijit._masterTT._onDeck=null;
//    }
//			
//    // cell tooltip
//    dojo.connect(grid, "onCellMouseOver", showTooltip);
//    dojo.connect(grid, "onCellMouseOut", hideTooltip);
//    // header cell tooltip
//    dojo.connect(grid, "onHeaderCellMouseOver", showTooltip);
//    dojo.connect(grid, "onHeaderCellMouseOut", hideTooltip);
//
//    // grid menu
//    window["gridMenu"] = dijit.byId("gridMenu");
//    gridMenu.bindDomNode(grid.domNode);
//    // prevent grid methods from killing the context menu event by implementing our own handler
//    grid.onCellContextMenu = function(e) {
//            cellNode = e.cellNode;
//    };
//    grid.onHeaderContextMenu = function(e) {
//        cellNode = e.cellNode;
//    };    
//    
//    reportCell = function() {
//        if(cellNode){
//            alert("Cell contents:  " + cellNode.innerHTML);
//            cellNode = null;
//        }
//    }
//    
//    gridTooltipEnabled = true;
//    toggleTooltip = function (button){
//        gridTooltipEnabled = !gridTooltipEnabled;
//        button.value = gridTooltipEnabled ? "Disable Grid Tooltip" : "Enable Grid Tooltip";
//    }
//    
//    gridMenuEnabled = true;
//    toggleMenu = function (button){
//        gridMenuEnabled = !gridMenuEnabled;
//        button.value = gridMenuEnabled ? "Disable Grid Menu" : "Enable Grid Menu";
//        gridMenu[gridMenuEnabled ? "bindDomNode" : "unBindDomNode"](grid.domNode);
//    }
//
//}
//
