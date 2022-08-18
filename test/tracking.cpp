#include <catch.hpp>

#include "databases.hpp"
#include "io.hpp"
#include "attitude-utils.hpp"

using namespace lost;

static unsigned char *BuildTrackingDatabase(const Catalog &catalog, long *length) {

    long dummyLength;
    if (length == NULL) length = &dummyLength;

    *length = SerializeLengthTrackingCatalog(catalog);
    unsigned char *result = new unsigned char[*length];
    SerializeTrackingCatalog(catalog, result);
    return result;
}

std::vector<int16_t> NaiveQuery(const Catalog catalog, const Vec3 point, float radius, float threshold) {
    std::vector<int16_t> correct_query_ind;
    for (int i = 0; i < (int)catalog.size(); i++) {
        CatalogStar cstar = catalog[i];
        Vec3 diff = cstar.spatial - point;
        if (diff.Magnitude() <= radius + threshold) {
            correct_query_ind.push_back(i);
        }
    }
    return correct_query_ind;
}

TEST_CASE("Tracking mode database", "[tracking]") {
    long length;
    Catalog &catalog = CatalogRead();
    unsigned char *dbBytes = BuildTrackingDatabase(catalog, &length);
    REQUIRE(length < 999999);
    TrackingSortedDatabase db(dbBytes);

    SECTION("small uncertainty") {

        for (int i = 0; i < (int)catalog.size(); i++) {
            Vec3 point = catalog[i].spatial;
            float radius = 0.001;
            float threshold = 0.001;

            std::vector<int16_t> query_ind = db.QueryNearestStars(catalog, point, radius, threshold);
            std::vector<int16_t> correct_query_ind = NaiveQuery(catalog, point, radius, threshold);

            std::sort(query_ind.begin(), query_ind.end());
            std::sort(correct_query_ind.begin(), correct_query_ind.end());
            CHECK(query_ind.size() == correct_query_ind.size());
            CHECK(query_ind == correct_query_ind);
        }
        
    }

    SECTION("medium uncertainty") {

        for (int i = 0; i < (int)catalog.size(); i++) {
            Vec3 point = catalog[i].spatial;
            float radius = 0.05;
            float threshold = 0.001;

            std::vector<int16_t> query_ind = db.QueryNearestStars(catalog, point, radius, threshold);
            std::vector<int16_t> correct_query_ind = NaiveQuery(catalog, point, radius, threshold);

            std::sort(query_ind.begin(), query_ind.end());
            std::sort(correct_query_ind.begin(), correct_query_ind.end());
            CHECK(query_ind.size() == correct_query_ind.size());
            CHECK(query_ind == correct_query_ind);
        }
    }

    SECTION("large uncertainty") {
        for (int i = 0; i < (int)catalog.size(); i++) {
            Vec3 point = catalog[i].spatial;
            float radius = 1;
            float threshold = 0.001;

            std::vector<int16_t> query_ind = db.QueryNearestStars(catalog, point, radius, threshold);
            std::vector<int16_t> correct_query_ind = NaiveQuery(catalog, point, radius, threshold);

            std::sort(query_ind.begin(), query_ind.end());
            std::sort(correct_query_ind.begin(), correct_query_ind.end());
            CHECK(query_ind.size() == correct_query_ind.size());
            CHECK(query_ind.size() != 0);
            CHECK(query_ind == correct_query_ind);
        }
    }
}

