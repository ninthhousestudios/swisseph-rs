/// Verify that every audited public-API type is importable from the crate root.
#[test]
fn imports_compile() {
    use swisseph::{
        // houses
        AscMc,
        // types — core
        AsteroidId,
        // types — astro models
        AstroModels,
        // azalt
        AzAltDir,
        BiasModel,
        Body,
        // flags
        CalcFlags,
        // context
        CalcResult,
        CalendarType,
        DegreeParts,
        DeltaT,
        DeltaTModel,
        EclipseFlags,
        // eclipse
        EclipseHow,
        EclipseWhere,
        Ephemeris,
        // config
        EphemerisConfig,
        EphemerisSource,
        Epsilon,
        // error
        Error,
        FictitiousBody,
        FictitiousId,
        FrameTransform,
        // heliacal
        HeliacalAngleResult,
        HeliacalEvent,
        HeliacalEventType,
        HeliacalFlags,
        HeliacalPheno,
        HorDir,
        HouseResult,
        HouseSystem,
        JdTt,
        JdUt1,
        JplHorMode,
        JplHoraMode,
        LunarEclipseGlobal,
        LunarEclipseHow,
        LunarEclipseLocal,
        // crossings
        MoonCrossing,
        // nodaps
        NodApsMethod,
        NodesApsides,
        Nutation,
        NutationModel,
        OccultGlobal,
        OccultLocal,
        // orbit
        OrbitalElements,
        // phenomena
        Phenomena,
        PlanetMoonId,
        PrecessionDirection,
        PrecessionModel,
        RefracDir,
        // crate-level alias
        Result,
        RiseSetFlags,
        // riseset
        RiseSetResult,
        SiderealMode,
        SiderealTimeModel,
        SolarEclipseGlobal,
        SolarEclipseLocal,
        SplitDegFlags,
        // stars
        Star,
        StarCatalog,
        TopoPosition,
        UtcComponents,
        UtcToJd,
        VisLimFlags,
        VisLimitResult,
    };

    // Touch each type to suppress unused-import warnings.
    let _ = std::any::type_name::<EphemerisConfig>();
    let _ = std::any::type_name::<TopoPosition>();
    let _ = std::any::type_name::<CalcResult>();
    let _ = std::any::type_name::<Ephemeris>();
    let _ = std::any::type_name::<Error>();
    let _ = std::any::type_name::<CalcFlags>();
    let _ = std::any::type_name::<EclipseFlags>();
    let _ = std::any::type_name::<HeliacalFlags>();
    let _ = std::any::type_name::<RiseSetFlags>();
    let _ = std::any::type_name::<SplitDegFlags>();
    let _ = std::any::type_name::<VisLimFlags>();
    let _ = std::any::type_name::<AsteroidId>();
    let _ = std::any::type_name::<Body>();
    let _ = std::any::type_name::<CalendarType>();
    let _ = std::any::type_name::<DegreeParts>();
    let _ = std::any::type_name::<EphemerisSource>();
    let _ = std::any::type_name::<Epsilon>();
    let _ = std::any::type_name::<FictitiousBody>();
    let _ = std::any::type_name::<FictitiousId>();
    let _ = std::any::type_name::<FrameTransform>();
    let _ = std::any::type_name::<HouseSystem>();
    let _ = std::any::type_name::<JdTt>();
    let _ = std::any::type_name::<JdUt1>();
    let _ = std::any::type_name::<Nutation>();
    let _ = std::any::type_name::<PlanetMoonId>();
    let _ = std::any::type_name::<PrecessionDirection>();
    let _ = std::any::type_name::<SiderealMode>();
    let _ = std::any::type_name::<UtcComponents>();
    let _ = std::any::type_name::<UtcToJd>();
    let _ = std::any::type_name::<AstroModels>();
    let _ = std::any::type_name::<BiasModel>();
    let _ = std::any::type_name::<DeltaTModel>();
    let _ = std::any::type_name::<JplHorMode>();
    let _ = std::any::type_name::<JplHoraMode>();
    let _ = std::any::type_name::<NutationModel>();
    let _ = std::any::type_name::<PrecessionModel>();
    let _ = std::any::type_name::<SiderealTimeModel>();
    let _ = std::any::type_name::<AzAltDir>();
    let _ = std::any::type_name::<HorDir>();
    let _ = std::any::type_name::<RefracDir>();
    let _ = std::any::type_name::<HeliacalAngleResult>();
    let _ = std::any::type_name::<HeliacalEvent>();
    let _ = std::any::type_name::<HeliacalEventType>();
    let _ = std::any::type_name::<HeliacalPheno>();
    let _ = std::any::type_name::<VisLimitResult>();
    let _ = std::any::type_name::<EclipseHow>();
    let _ = std::any::type_name::<EclipseWhere>();
    let _ = std::any::type_name::<LunarEclipseGlobal>();
    let _ = std::any::type_name::<LunarEclipseHow>();
    let _ = std::any::type_name::<LunarEclipseLocal>();
    let _ = std::any::type_name::<OccultGlobal>();
    let _ = std::any::type_name::<OccultLocal>();
    let _ = std::any::type_name::<SolarEclipseGlobal>();
    let _ = std::any::type_name::<SolarEclipseLocal>();
    let _ = std::any::type_name::<AscMc>();
    let _ = std::any::type_name::<HouseResult>();
    let _ = std::any::type_name::<NodApsMethod>();
    let _ = std::any::type_name::<NodesApsides>();
    let _ = std::any::type_name::<OrbitalElements>();
    let _ = std::any::type_name::<Phenomena>();
    let _ = std::any::type_name::<RiseSetResult>();
    let _ = std::any::type_name::<MoonCrossing>();
    let _ = std::any::type_name::<Star>();
    let _ = std::any::type_name::<StarCatalog>();
    let _ = std::any::type_name::<Result<()>>();

    // DeltaT is a trait — verify it's usable as a bound.
    fn _assert_delta_t<T: DeltaT>(_t: &T) {}
}
