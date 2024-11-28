function orphan(model::ModelSource{ID, Dom}) where ID where Dom
    P = typeof(model.parameters)
    DummyModel{ID, Dom, P}(model.parameters)
end

function orphan(model::Dolo.YModel)
    Dolo.YModel(
        Dolo.name(model),
        model.states,
        model.controls,
        model.exogenous,
        model.calibration,
        orphan(model.source)
    )
end

function orphan(model::Dolo.DYModel)
    Dolo.DYModel(
        orphan(model.model),
        model.grid,
        model.dproc
    )
end

